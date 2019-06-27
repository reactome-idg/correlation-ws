package org.reactome.idg.loader;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Calculates a matrix of N x N (where N is number of genes) for a specific tissue.
 * @author sshorser
 *
 */
public class GenePairCorrelationMatrixCalculator extends CorrelationCalculator
{
	private static final Logger logger = LogManager.getLogger();
	
	/**
	 * Creates a calculator.
	 * @param t - the path to the file with the tissue samples.
	 * @param loader - the data loader, associated with a specific HDF file.
	 */
	public GenePairCorrelationMatrixCalculator(String t, Archs4ExpressionDataLoader loader)
	{
		super(t, loader);
	}
	
	/**
	 * Calculate the correlation between all genes for a specific tissue.
	 * @return a triangular matrix containing correlation values. Only the main diagonal and values above will be populated (the whole thing does NOT need to be populated since the lower half is just the
	 * mirror image of the upper half). The matrix will be N x N, where N = number of genes in HDF file.
	 * @throws IOException
	 */
	public RealMatrix calculateCorrelation() throws IOException
	{
		// We're assuming that "tissues" are names of files with tissue-specific sample IDs.
		// outer index is sample, inner index is gene.
		int[][] sampleValues = this.dataLoader.getExpressionValuesforTissue(Paths.get(this.tissue));
		int numberOfGenes = this.dataLoader.getGeneIndices().size();
//		double corMatrix[][] = new double[numberOfGenes][numberOfGenes];
		RealMatrix corMatrix = new Array2DRowRealMatrix(numberOfGenes, numberOfGenes);
		ForkJoinPool pool = new ForkJoinPool();
		List<Callable<Double>> workers = new ArrayList<>();
		AtomicInteger calculationsPerformed = new AtomicInteger(0);
		// need to compare each gene to each other gene.
		for (int geneIndex = 0; geneIndex < numberOfGenes - 1; geneIndex++)
		{
			int numberOfSamples = sampleValues.length;
			// Build a list of sample values for a gene, for all the samples associated with the given tissue.
			final double[] geneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, geneIndex, numberOfSamples);
			// Now, for all OTHER genes...
			for (int otherGeneIndex = geneIndex; otherGeneIndex < numberOfGenes - 1; otherGeneIndex ++)
			{
				final double[] otherGeneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, otherGeneIndex, numberOfSamples);
				// now that we have a list of values for a gene, and for another gene, we can calculate the correlation:
				final int gIndx = geneIndex;
				final int gOtherIndx = otherGeneIndex;
				// Create a worker than will calculate the correlation for two lists of sample data.
				Callable<Double> worker = new Callable<Double>()
				{
					@Override
					public Double call() throws Exception
					{
						PearsonsCorrelation cor = new PearsonsCorrelation();
						double correlationValue = cor.correlation(geneSamples, otherGeneSamples);
						calculationsPerformed.getAndIncrement();
//						synchronized (corMatrix)
						{
							corMatrix.setEntry(gIndx, gOtherIndx, correlationValue);
						}
						// data gets saved into a shared array, so the return value won't really go anywhere.
						return correlationValue;
					}
				};
				workers.add(worker);
				// When there are enough workers waiting, invoke them all. 10,000 seems to be the most optimal group size,
				// on my workstation (8 cores x @~3GHz, 64 GB RAM)
				if (workers.size() % 5000 == 0)
				{
//					logger.info("# workers: " + workers.size());
					pool.invokeAll(workers);
					workers.clear();
				}
			}

			if (geneIndex % 1000 == 0)
			{
				logger.info(calculationsPerformed.get() + " calculations performed. Current Gene Index: " + geneIndex);
			}
		}
		if (workers.size() > 0)
		{
			logger.info("Invoking {} remaining workers...", workers.size());
			pool.invokeAll(workers);
			workers.clear();
		}
		
		return corMatrix;
	}
}
