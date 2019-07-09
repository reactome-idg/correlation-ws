package org.reactome.idg.verification;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.DataRepository;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;
import org.reactome.idg.loader.Archs4ExpressionDataLoaderFactory;
import org.reactome.idg.loader.GenePairCorrelationCalculator;

public class CorrelationVerification
{
	private static final Logger logger = LogManager.getLogger();
	
	public static void main(String[] args) throws IOException
	{
		String pathToHDFFile = args[0];
		int numPairs = Integer.parseInt(args[1]);
		Map<String, Double> calculatedCorrelations = Collections.synchronizedMap(new HashMap<>(numPairs));
		// First, we need to read the HDF file.
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(pathToHDFFile);
		loader.loadMetaData();
		StringBuffer buf = new StringBuffer();
		buf.append("Gene Pair").append("\t").append("Old Correlation (extracted from correlation file)").append("\t").append("New Correlation (calculated from expression data in HDF)\n");
		logger.info("Number of correlations to calculate: {}", numPairs);
		String pathToCorrelationFile = args[2];
		Map<String, Integer> geneNamesToIndices = new HashMap<>();
		// read the file. Start by reading the header to get the indices.
		try(FileInputStream fis = new FileInputStream(pathToCorrelationFile);
				Scanner scanner = new Scanner(fis);)
		{
			if (scanner.hasNext())
			{
				String geneSymbolHeader = scanner.nextLine();
				String[] parts = geneSymbolHeader.split(",");
				String[] allGeneSymbols = new String[parts.length-1];
				for (int i = 1; i < parts.length; i++)
				{
					allGeneSymbols[i - 1] = parts[i].replaceAll("\"", "");
					geneNamesToIndices.put(allGeneSymbols[i - 1], i);
				}

				logger.info("{} gene symbols in the header.", allGeneSymbols.length);
			}
		}
		
		Random random = new Random(); // should I bother seeding this with the current time?
		// Need to know how large the random numbers can be.
		int maxNum = loader.getGeneIndices().keySet().size();
		// get a bunch of random numbers
//		for (int i : random.ints(numPairs, 0, maxNum).toArray())
		random.ints(numPairs, 0, maxNum).parallel().forEach( i -> 
		{
			int geneIndex1 = i;
			String geneSymbol1 = loader.getGeneIndicesToNames().get(geneIndex1);
			while (!geneNamesToIndices.containsKey(geneSymbol1))
			{
				geneIndex1 = random.nextInt(maxNum);
				geneSymbol1 = loader.getGeneIndicesToNames().get(geneIndex1);
			}
			String geneSymbol2 = "";
			int otherIndex = random.nextInt(maxNum);
			// Make sure we don't select the same two genes, or the correlation will just be 1.0, and that's not very interesting.
			// Actually, maybe it should be allowed. If two random genes ARE the same, then the output should be 1.0 - if it's NOT 
			// then that would indicate something else has gone wrong.
			while (geneSymbol2 == geneSymbol1 || !geneNamesToIndices.containsKey(geneSymbol2))
			{
				otherIndex = random.nextInt(maxNum);
				geneSymbol2 = loader.getGeneIndicesToNames().get(otherIndex);
			}
			logger.info("Two genes selected: {} and {}", geneSymbol1, geneSymbol2);
			// Now that we have two genes, compute cross-tissue correlations.
			GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator(geneSymbol1, geneSymbol2, loader);
			double correlation;
			try
			{
//				synchronized(calculatedCorrelations)
				{
					// Need to do the calculation in a synchronized block, because reading from HDF cannot be safely done in a multi-threaded way.
					correlation = calculator.calculateGenePairCorrelation();
					logger.info("Correlation between {} and {} is: {}", geneSymbol1, geneSymbol2, correlation);
					calculatedCorrelations.put(DataRepository.generateKey(geneSymbol1, geneSymbol2), correlation);
				}
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		});
		logger.info("Now getting correlations values from file.");
		Map<String,Double> correlationsFromFile = Collections.synchronizedMap(new HashMap<>(numPairs));
		LinkedBlockingQueue<String> lines = new LinkedBlockingQueue<>();
		
		ExecutorService execService = Executors.newCachedThreadPool();
		for (int i = 0; i < ForkJoinPool.getCommonPoolParallelism(); i++)
		{
			Runnable lineProcessor = new Runnable()
			{
				@Override
				public void run()
				{
					while (correlationsFromFile.size() != calculatedCorrelations.size())
					{
//						if (lines.poll() != null)
						{
							String line;
							try
							{
								// the take() method is blocking, so this Runnable *will* throw an InterruptedException when the executor service is shut down.
								line = lines.take();
								String gene1 = calculatedCorrelations.keySet().parallelStream().map( geneKey -> geneKey.split("\\|")[0] ).filter( gene -> line.startsWith("\""+gene+"\"") ).findFirst().orElse(null) ;
								if (gene1 != null)
								{
									String geneKey = calculatedCorrelations.keySet().parallelStream().filter(gk -> gk.startsWith(gene1 + "|")).findFirst().get();
									String[] parts = geneKey.split("\\|");
									String gene2 = parts[1];
									int geneIndex2 = geneNamesToIndices.get(gene2);
									// If the current line is for gene 1, get the correlation from the line.
									if (line.startsWith("\""+gene1+"\"") )
									{
										// now we found the line we need, get the correlation value.
										String[] lineParts;
										lineParts = line.split(",");
										double correlationFromFile = Double.parseDouble(lineParts[geneIndex2]);
										logger.info("Correlation (according to the FILE) between {} and {} is: {}", gene1, gene2, correlationFromFile);
										correlationsFromFile.put(geneKey, correlationFromFile);
										logger.info("Delta: {}", calculatedCorrelations.get(geneKey) - correlationsFromFile.get(geneKey) );
										buf.append(geneKey.replace("|", ", ")).append("\t").append(correlationsFromFile.get(geneKey)).append("\t").append(calculatedCorrelations.get(geneKey)).append("\n");
									}
								}
							}
							catch (InterruptedException e)
							{
								logger.info("Interrupted.");
							}
						}
					}
				}
			};
			execService.submit(lineProcessor);
		}
		// now that we have the indices, we need to read into the file.
		// need to find a way to read to a specific line number in a file, using Java.
		try(FileInputStream fis = new FileInputStream(pathToCorrelationFile);
				Scanner scanner = new Scanner(fis);)
		{
			while (scanner.hasNext() && correlationsFromFile.size() != calculatedCorrelations.size())
			{
				String line = scanner.nextLine();
				lines.add(line);
//				String gene1 = calculatedCorrelations.keySet().parallelStream().map( geneKey -> geneKey.split("\\|")[0] ).filter( gene -> line.startsWith("\""+gene+"\"") ).findFirst().orElse(null) ;
//				if (gene1 != null)
//				{
//					String geneKey = calculatedCorrelations.keySet().parallelStream().filter(gk -> gk.startsWith(gene1 + "|")).findFirst().get();
//					String[] parts = geneKey.split("\\|");
//					String gene2 = parts[1];
//					int geneIndex2 = geneNamesToIndices.get(gene2);
//					// If the current line is for gene 1, get the correlation from the line.
//					if (line.startsWith("\""+gene1+"\"") )
//					{
//						// now we found the line we need, get the correlation value.
//						String[] lineParts;
//						lineParts = line.split(",");
//						double correlationFromFile = Double.parseDouble(lineParts[geneIndex2]);
//						logger.info("Correlation (according to the FILE) between {} and {} is: {}", gene1, gene2, correlationFromFile);
//						correlationsFromFile.put(geneKey, correlationFromFile);
//						logger.info("Delta: {}", calculatedCorrelations.get(geneKey) - correlationsFromFile.get(geneKey) );
//						buf.append(geneKey.replace("|", ", ")).append("\t").append(correlationsFromFile.get(geneKey)).append("\t").append(calculatedCorrelations.get(geneKey)).append("\n");
//						if (correlationsFromFile.size() == calculatedCorrelations.size())
//						{
//							done = true;
//						}
//					}
//				}
			}
			execService.shutdownNow();
		}
		
		Files.write(Paths.get("./correlations_comparison.tsv"), buf.toString().getBytes(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.WRITE);
	}
}
