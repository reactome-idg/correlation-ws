package org.reactome.idg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.FileWriter;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.temporal.Temporal;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;
import org.reactome.idg.loader.Archs4ExpressionDataLoaderFactory;
import org.reactome.idg.loader.GenePairCorrelationCalculator;
import org.reactome.idg.loader.GenePairCorrelationMatrixCalculator;



@SuppressWarnings("static-method")
public class TestCorrelationCalculator
{

	private static final String PATH_TO_HDF = "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5";
	@Test
	public void testCorrelationCalculationsIT() throws IOException
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		// It takes about 20 mintues to calculate all gene pair-wise correlations for the heart.txt samples, 141 samples for all genes (~35k) in the H5 file.
		GenePairCorrelationMatrixCalculator heartCorMatrixCalculator = new GenePairCorrelationMatrixCalculator("src/test/resources/heart.txt", loader);
		
		RealMatrix m = heartCorMatrixCalculator.calculateCorrelation();
		for (int i = 0; i < Math.min(m.getColumnDimension(), 20); i++)
		{
			for (int j = 0; j < Math.min(m.getRowDimension(), 20); j++)
			{
				System.out.print(m.getEntry(i, j) + "\t");
			}
			System.out.println();
		}
	}
	
	@Test
	public void testGenePairCorrelationsIT() throws IOException
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator("src/test/resources/heart.txt", "A1BG","A1BG", loader);
		double d = calculator.calculateGenePairCorrelation();
		System.out.println(d);
		
		assertEquals(d, 1.0, 0);
		
		d = calculator.calculateGenePairCorrelation("A1BG","A1CF");
		System.out.println(d);
		
		try
		{
			d = calculator.calculateGenePairCorrelation("BLAHBLAHBLAH","A1CF");
			System.out.println(d);
		}
		catch (IllegalArgumentException e)
		{
			e.printStackTrace();
			assertTrue(e.getMessage().contains("Gene BLAHBLAHBLAH is not recognized in the HDF file."));
		}
		
		try
		{
			d = calculator.calculateGenePairCorrelation("A1CF", "BLAH6666666");
			System.out.println(d);
		}
		catch (IllegalArgumentException e)
		{
			e.printStackTrace();
			assertTrue(e.getMessage().contains("Gene BLAH6666666 is not recognized in the HDF file."));
		}
	}
	

	/*
	 * JVM args used:
	 * -ea -Xms20G -XX:+UseParallelGC
	 */
	@Test
	public void testGenePairCalculatorAllGenesIT() throws InterruptedException, IOException, OutOfRangeException, ExecutionException
	{
		// create a gene-pair list and then run GenePairCorrelationCalculator for each pair. See if it's an faster than the other method.
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		int numGenes = loader.getGeneIndices().keySet().size();
		int numPairs = ((1 + loader.getGeneIndices().keySet().size()) * loader.getGeneIndices().keySet().size())/2;
		System.out.println(numPairs + " possible gene-pairs for " + loader.getGeneIndices().keySet().size() + " genes.");
		int i = 0;
		// RealMatrix is a little faster to use in this context than double[][] 
		RealMatrix cor = new Array2DRowRealMatrix(numGenes, numGenes);
//		try(FileWriter fileWriter = new FileWriter("/media/sshorser/data/correlations_test", true);)
		{
			System.out.println("getCommonPoolParallelism: " + ForkJoinPool.getCommonPoolParallelism());
			int queueCapacity = 10000000 ; //ForkJoinPool.getCommonPoolParallelism() * 1000;
			System.out.println("Queue capacity: " + queueCapacity);
			BlockingQueue<Runnable> jobQueue = new ArrayBlockingQueue<>(queueCapacity);
			ThreadPoolExecutor execService = new ThreadPoolExecutor(ForkJoinPool.getCommonPoolParallelism() / 2, ForkJoinPool.getCommonPoolParallelism(), 10, TimeUnit.SECONDS, jobQueue);
			
			long totalDrainTime = 0;
			int numDrains = 0;
			System.out.println("building gene-pair list...");
//			GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator("src/test/resources/heart.txt", loader);
			for (int outer = 0; outer < numGenes; outer ++)
			{
				for (int inner = outer; inner < numGenes; inner++)
				{
					final int indx1 = outer;
					final int indx2 = inner;
					final String g1 = loader.getGeneIndicesToNames().get(indx1);
					final String g2 = loader.getGeneIndicesToNames().get(indx2);
					Callable<Double> calcWorker = new Callable<Double>()
					{
						@Override
						public Double call() throws Exception
						{
							try
							{
								GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator("src/test/resources/heart.txt", loader);
								double d = calculator.calculateGenePairCorrelation(g1, g2);
//								synchronized(fileWriter)
								{
									cor.setEntry(indx1, indx2, d);
//									fileWriter.write(g1+"|"+g2+"|"+d+"\n");
								}
								return d;
							}
							catch (IOException e)
							{
								e.printStackTrace();
								throw new RuntimeException(e);
							}
						}
					};
					// sleep when capacity gets low. Should prevent GC thrashing...
//					while (execService.getQueue().remainingCapacity() < 100)
//					{
//						Thread.sleep(1000);
//					}
					
//					if (execService.getQueue().size() % 100 == 0)
//					{
//						System.out.println("queue size: " + execService.getQueue().size());
//					}
					int remainingCap = execService.getQueue().remainingCapacity();
					int drainThreshold = Math.min(queueCapacity / 10, 20);
					
					if (remainingCap < drainThreshold)
					{
//						System.out.println("Remaining capacity: " + remainingCap + " drain threshold: " + drainThreshold);
						LocalDateTime drainStartTime = LocalDateTime.now();
						// sleep until the queue has room in it.
						while (execService.getQueue().size() > 0)
						{
							Thread.sleep(1);
//							System.out.println("queue size: " + execService.getQueue().size());
						}
						LocalDateTime drainEndTime = LocalDateTime.now();
						numDrains++;
						totalDrainTime += Duration.between(drainStartTime, drainEndTime).toMillis();
						String avgDrainTime = Duration.ofMillis(totalDrainTime / numDrains).toString();
//						System.out.println("Drain time: " + Duration.between(drainStartTime, drainEndTime) + " Avg drain time: " + avgDrainTime);
					}
					execService.submit(calcWorker);
//					cor.setEntry(indx1, indx2, execService.submit(calcWorker).get());

//					calculators.forEach(c -> c.run());
//					genePairQueue.notify();
					i++;
//					if (i % 100000 == 0)
//					{
//						Random r = new Random();
//						System.out.println("quick rest "+i+" "+cor[r.nextInt(outer)][r.nextInt(inner)]);
//						Thread.currentThread().sleep(1000);
//					}
//					if (i % 1000000 == 0)
//					{
//
//						float percentComplete = (float) (i * 1.0 / numPairs * 100.0);
//						System.out.println(LocalDateTime.now().toString() + "\t\t" + i/1000000 + " M ("+percentComplete+"% of all pairs)\t" + " queue size: " + execService.getQueue().size());
//					}
					// if the genePairQ gets too big, rest a moment so that the workers can catch up. Otherwise, the system runs slower and slower....
/*
					if (genePairQueue.size() > 1000)
					{
						LocalDateTime beforeRestEnd = LocalDateTime.now();
//						System.out.println("Time since last rest " + Duration.between(beforeRestStart, beforeRestEnd).toString());
						
						LocalDateTime start = LocalDateTime.now();
//						System.out.println("Letting GenePair queue drain...");
						while (genePairQueue.size() > 0)
						{
							Thread.sleep(50);
						}
						LocalDateTime end = LocalDateTime.now();
//						System.out.println("Elapsed time to drain queue: " + Duration.between(start, end).toString());
						beforeRestStart = LocalDateTime.now();
					}
*/
				}
				if (outer % 1000 == 0)
				{
					System.out.println(LocalDateTime.now().toString() + "\tGene index: " + outer);
				}
				
			}
//			calculators.forEach(c->c.end());
//			writers.forEach(w->w.end());
//			calcExecService.shutdown();
//			writeExecService.shutdown();
			execService.shutdown();
		}
	}
}
