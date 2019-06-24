package org.reactome.idg.loader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public abstract class CorrelationCalculator
{
	protected String tissue;
	protected Archs4ExpressionDataLoader dataLoader ;

	public CorrelationCalculator(String t, Archs4ExpressionDataLoader loader)
	{
		this.tissue = t;
		this.dataLoader = loader;
	}
	
	protected static double[] getSampleValuesForGene(int[][] sampleValues, int geneIndex, int numberOfSamples )
	{
		double[] geneSamples = new double[numberOfSamples];
		for (int sampleIndex = 0; sampleIndex < numberOfSamples; sampleIndex++)
		{
			geneSamples[sampleIndex] = sampleValues[sampleIndex][geneIndex];
		}
		return geneSamples;
	}
	
}
