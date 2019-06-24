package org.reactome.idg.loader;

/**
 * @author sshorser
 *
 */
public abstract class CorrelationCalculator
{
	protected String tissue;
	protected Archs4ExpressionDataLoader dataLoader ;

	/**
	 * Create a new calculator for a specific tissue (whose sample IDs come from a file), using a specific Archs4 data loader.
	 * @param t - the path to the tissue sample file.
	 * @param loader - the data loader, associated with a specific HDF file.
	 */
	public CorrelationCalculator(String t, Archs4ExpressionDataLoader loader)
	{
		this.tissue = t;
		this.dataLoader = loader;
	}

	/**
	 * Gets the sample values for a specific gene, from a set of sample values.
	 * @param sampleValues - the list of ALL sample values.
	 * @param geneIndex - the index in the array of the gene.
	 * @param numberOfSamples - the total number of samples in the result.
	 * @return
	 */
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
