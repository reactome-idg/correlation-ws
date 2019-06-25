package org.reactome.idg.loader;

import java.io.IOException;
import java.nio.file.Paths;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GenePairCorrelationCalculator extends CorrelationCalculator
{

	private String gene1;
	private String gene2;
	private static final Logger logger = LogManager.getLogger();
	private String currentTissue;
	private static int[][] cachedExprValues;
	
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
		this.verifyGenes();
		return this.calculateGenePairCorrelation();
	}
	
	/**
	 * Calculate the correlation between the genes that were passed to the constructor.
	 * @return the correlation value.
	 * @throws IOException
	 */
	public double calculateGenePairCorrelation() throws IOException
	{
		this.verifyGenes();
		int[][] sampleValues = new int[1][1];
		if (!this.tissue.equals(currentTissue))
		{
			currentTissue = this.tissue;
			sampleValues = this.dataLoader.getExpressionValuesforTissue(Paths.get(this.tissue));
			cachedExprValues = sampleValues;
		}
		else
		{
			if (cachedExprValues != null)
			{
				sampleValues = cachedExprValues;
			}
		}
		int numberOfSamples = sampleValues.length;

		
		int geneIndex = this.dataLoader.getGeneIndices().get(this.gene1);
		int otherGeneIndex = this.dataLoader.getGeneIndices().get(this.gene2);
		
		final double[] geneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, geneIndex, numberOfSamples);
		final double[] otherGeneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, otherGeneIndex, numberOfSamples);
		
		
		
		PearsonsCorrelation cor = new PearsonsCorrelation();
		double correlationValue = cor.correlation(geneSamples, otherGeneSamples);
		
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
