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

public class GenePairCorrelationCalculator
{
	private static final Logger logger = LogManager.getLogger();
//	private String gene1;
//	private String gene2;
	private String tissue;
	
	private double[] gene1SampleValues;
	private double[] gene2SampleValues;
	
	public GenePairCorrelationCalculator(String t, String hdfExpressionFile)
	{
//		this.gene1 = g1;
//		this.gene2 = g2;
		this.tissue = t;
		Archs4ExpressionDataLoader.setHdfExpressionFileAndLoadCounts(hdfExpressionFile);
		Archs4ExpressionDataLoader.loadMetaData();
	}
	
	public double[][] calculateCorrelation() throws IOException
	{
		
		
		// We're assuming that "tissues" are names of files with tissue-specific sample IDs.
		// outer index is sample, inner index is gene.
		int[][] sampleValues = Archs4ExpressionDataLoader.getExpressionValuesforTissue(Paths.get(this.tissue));
		int numberOfGenes = Archs4ExpressionDataLoader.getGeneIndices().size();
		double corMatrix[][] = new double[numberOfGenes][numberOfGenes];
		ForkJoinPool pool = new ForkJoinPool();
		List<Callable<Double>> workers = new ArrayList<>();
		AtomicInteger calculationsPerformed = new AtomicInteger(0);
		// need to compare each gene to each other gene.
		for (int geneIndex = 0; geneIndex < numberOfGenes - 1; geneIndex++)
		{
			int numberOfSamples = sampleValues.length;
			double[] geneSamples = new double[numberOfSamples];
			for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++)
			{
				geneSamples[sampleIndex] = sampleValues[sampleIndex][geneIndex];
			}
			for (int otherGeneIndex = geneIndex; otherGeneIndex < numberOfGenes - 1; otherGeneIndex ++)
			{
				double[] otherGeneSamples = new double[numberOfSamples];
				for (int otherSampleIndex = 0; otherSampleIndex < numberOfSamples; otherSampleIndex++)
				{
					otherGeneSamples[otherSampleIndex] = sampleValues[otherSampleIndex][otherGeneIndex];
				}
				// now that we have a list of values for a gene, and for another gene, we can calculate the correlation:
				String gene1, gene2;
				gene1 = Archs4ExpressionDataLoader.getGeneIndicesToNames().get(geneIndex);
				gene2 = Archs4ExpressionDataLoader.getGeneIndicesToNames().get(otherGeneIndex);

				final int gIndx = geneIndex;
				final int gOtherIndx = otherGeneIndex;
				Callable<Double> worker = new Callable<Double>()
				{
					@Override
					public Double call() throws Exception
					{
						PearsonsCorrelation cor = new PearsonsCorrelation();
						double correlationValue = cor.correlation(geneSamples, otherGeneSamples);
						calculationsPerformed.getAndIncrement();
						synchronized (corMatrix)
						{
							corMatrix[gIndx][gOtherIndx] = correlationValue;
						}
						
						return correlationValue;
					}
				};
				workers.add(worker);
				if (workers.size() % 10000 == 0)
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
