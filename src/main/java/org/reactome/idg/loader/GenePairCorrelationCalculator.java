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
	
	public GenePairCorrelationCalculator(String t, String gene1, String gene2, Archs4ExpressionDataLoader loader)
	{
		super(t, loader);
		this.gene1 = gene1;
		this.gene2 = gene2;
	}
	
	public double calculateGenePairCorrelation(String g1, String g2) throws IOException
	{
		this.gene1 = g1;
		this.gene2 = g2;
		return this.calculateGenePairCorrelation();
	}
	
	public double calculateGenePairCorrelation() throws IOException
	{
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

		if (!this.dataLoader.getGeneIndices().containsKey(gene1))
		{
			throw new IllegalArgumentException("Gene "+gene1+" is not recognized in the HDF file.");
		}
		if (!this.dataLoader.getGeneIndices().containsKey(gene2))
		{
			throw new IllegalArgumentException("Gene "+gene2+" is not recognized in the HDF file.");
		}
		
		int geneIndex = this.dataLoader.getGeneIndices().get(this.gene1);
		int otherGeneIndex = this.dataLoader.getGeneIndices().get(this.gene2);
		
		final double[] geneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, geneIndex, numberOfSamples);
		final double[] otherGeneSamples = CorrelationCalculator.getSampleValuesForGene(sampleValues, otherGeneIndex, numberOfSamples);
		
		
		
		PearsonsCorrelation cor = new PearsonsCorrelation();
		double correlationValue = cor.correlation(geneSamples, otherGeneSamples);
		
		return correlationValue;
	}
}
