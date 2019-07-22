package org.reactome.idg.loader;

import java.io.IOException;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.util.Arrays;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Can be used to calculate the correlation between two specific genes, for a specific tissue.
 * @author sshorser
 *
 */
public class GenePairCorrelationCalculator extends CorrelationCalculator
{

	private String gene1;
	private String gene2;
	private static final Logger logger = LogManager.getLogger();
	private String currentTissue;
	private static double[][] cachedExprValues;
	
	/**
	 * Creates a calculator without specify the genes.
	 * @param t
	 * @param loader
	 */
	public GenePairCorrelationCalculator(String t, Archs4ExpressionDataLoader loader)
	{
		super(t, loader);
	}
	
	/**
	 * Creates a new calculator that will calculate the correlation for a singel gene pair.
	 * @param t - the path to the tissue sample file.
	 * @param gene1 - the first gene.
	 * @param gene2 - the second gene.
	 * @param loader - an Archs4ExpressionDataLoader that is associated with a specific HDF file.
	 */
	public GenePairCorrelationCalculator(String t, String gene1, String gene2, Archs4ExpressionDataLoader loader)
	{
		super(t, loader);
		this.gene1 = gene1;
		this.gene2 = gene2;
		this.verifyGenes();
	}
	
	/**
	 * Creates a correlation calculator that will use ALL expression data available. The calculated value will NOT be tissue-specific.
	 * @param gene1 - the first gene.
	 * @param gene2 - the second gene.
	 * @param loader - a Loader for a specific HDF file.
	 */
	public GenePairCorrelationCalculator(String gene1, String gene2, Archs4ExpressionDataLoader loader)
	{
		super(null, loader);
		this.gene1 = gene1;
		this.gene2 = gene2;
		this.verifyGenes();
	}
	
	/**
	 * This will set gene1 and gene2 for this calculator and then calculate the correlation.
	 * @param g1
	 * @param g2
	 * @return the correlation value.
	 * @throws IOException
	 */
	public double calculateGenePairCorrelation(String g1, String g2) throws IOException
	{
		this.gene1 = g1;
		this.gene2 = g2;
		return this.calculateGenePairCorrelation();
	}
	
	/**
	 * Calculate the correlation between the genes that were passed to the constructor. If tissue is not specified, correlation will be calculated across ALL samples.
	 * @return the correlation value.
	 * @throws IOException
	 */
	public double calculateGenePairCorrelation() throws IOException
	{
		this.verifyGenes();
		
		if (this.tissue == null)
		{
			return calculateCrossTissueCorrelation();
		}
		else
		{
			double[][] sampleValues = new double[1][1];
			// if the tissue has changed, we'll need to load new values into the cache.
			if (!this.tissue.equals(currentTissue))
			{
				currentTissue = this.tissue;
				sampleValues = this.dataLoader.getExpressionValuesforTissue(Paths.get(this.tissue));	
				
				cachedExprValues = sampleValues;
			}
			else
			{
				// if there is a cache then use it.
				if (cachedExprValues != null)
				{
					sampleValues = cachedExprValues;
				}
				else 
				{
					// this code path probably isn't possible, since the current and previous tissues are already known to match, meaning the samples have probaby already been loaded.
					sampleValues = this.dataLoader.getExpressionValuesforTissue(Paths.get(this.tissue));
					cachedExprValues = sampleValues;
				}
			}
			int numberOfSamples = sampleValues.length;
			// get indices of the genes.
			int geneIndex = this.dataLoader.getGeneIndices().get(this.gene1);
			int otherGeneIndex = this.dataLoader.getGeneIndices().get(this.gene2);
			// get the sample values for the two genes.
			final double[] geneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, geneIndex, numberOfSamples);
			final double[] otherGeneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, otherGeneIndex, numberOfSamples);
			// calculate the Pearson's Correlation for the two sets of values, and return the correlation value.
			PearsonsCorrelation cor = new PearsonsCorrelation();
			double correlationValue = cor.correlation(geneSamples, otherGeneSamples);
			return correlationValue;
		}
	}
	
	public double calculateCrossTissueCorrelation()
	{
		double[] gene1ExpressionValues;
		double[] gene2ExpressionValues;
		gene1ExpressionValues = this.dataLoader.getAllExpressionValuesForGene(this.gene1);
		gene2ExpressionValues = this.dataLoader.getAllExpressionValuesForGene(this.gene2);
		
		PearsonsCorrelation cor = new PearsonsCorrelation();
		double correlationValue = cor.correlation(Arrays.stream(gene1ExpressionValues).toArray(),
												Arrays.stream(gene2ExpressionValues).toArray());
		
		return correlationValue;
	}

//	private static double[] applyJLTransform(int[] geneExpressionValues)
//	{
//		JLTransform transform1 = new JLTransform(1000);
//		List<DataPoint> dataPoints1 = new ArrayList<>();
//		Vec vec1 = Vec.zeros(geneExpressionValues.length);
//		for (int i = 0; i < geneExpressionValues.length; i++)
//		{
//			vec1.set(i, geneExpressionValues[i]);
//		}
//		dataPoints1.add(new DataPoint(vec1 ));
//		DataSet<?> ds1 = new SimpleDataSet(dataPoints1);
//		transform1.fit(ds1);
//		ds1.applyTransform(transform1);
//		List<DataPoint> postTransformDatapoints = ds1.getDataPoints();
//		Vec vec2 = postTransformDatapoints.get(0).getNumericalValues();
//		double[] transformedValues = new double[vec2.length()];
//		for (int i = 0; i < vec2.length(); i++)
//		{
//			transformedValues[i] = vec2.get(i);
//		}
//		return transformedValues;
//	}
//	
//	private double[] normalizeExpressionValues(int[] expressionValues, int geneIndex)
//	{
//		double[] normalized = new double[expressionValues.length];
//		
//		for (int i = 0; i < expressionValues.length; i++)
//		{
//			// Now we need to create normalized sample values, then only keep the one for the gene in question.
//			// Normalizing the WHOLE matrix will not fit in memory, so we need to do it this way.
//			int[][] expr = this.dataLoader.getExpressionValuesBySampleIndices(Arrays.asList(i), "");
//			double[] norm = StatUtils.normalize(Arrays.stream(expr[0]).asDoubleStream().toArray());
//			normalized[i] = norm[geneIndex];
//		}
//		
//		return normalized;
//	}
	
	public double calculateCrossTissueCorrelationFilteredByDate(LocalDateTime cutoff)
	{
		double[] gene1ExpressionValues;
		double[] gene2ExpressionValues;
		gene1ExpressionValues = this.dataLoader.getDateFilteredExpressionValuesForGene(this.gene1, cutoff);
		gene2ExpressionValues = this.dataLoader.getDateFilteredExpressionValuesForGene(this.gene2, cutoff);
//		int gene1Index = this.dataLoader.getGeneIndices().get(this.gene1);
//		double[] norm1 = normalizeExpressionValues(gene1ExpressionValues, gene1Index);
//		
//		int gene2Index = this.dataLoader.getGeneIndices().get(this.gene1);
//		double[] norm2 = normalizeExpressionValues(gene2ExpressionValues, gene2Index);
//		
		PearsonsCorrelation cor = new PearsonsCorrelation();
//		double correlationValue = cor.correlation(norm1, norm2);
		double correlationValue = cor.correlation(StatUtils.normalize(Arrays.stream(gene1ExpressionValues).toArray()),
												StatUtils.normalize(Arrays.stream(gene2ExpressionValues).toArray()));
		
		return correlationValue;
	}
	
	/**
	 * Verify that the genes are OK. If they are not known in the HDF file, then thrown an exception.
	 * @throws IllegalArgumentException If genes are not valid.
	 */
	private void verifyGenes()
	{
		if (!this.dataLoader.getGeneIndices().containsKey(gene1))
		{
			String message = "Gene "+gene1+" is not recognized in the HDF file.";
			logger.error(message);
			throw new IllegalArgumentException(message);
		}
		if (!this.dataLoader.getGeneIndices().containsKey(gene2))
		{
			String message = "Gene "+gene1+" is not recognized in the HDF file.";
			logger.error(message);
			throw new IllegalArgumentException(message);
		}
	}

}
