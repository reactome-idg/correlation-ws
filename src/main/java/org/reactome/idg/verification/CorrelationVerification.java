package org.reactome.idg.verification;

import java.io.IOException;
import java.util.Random;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;
import org.reactome.idg.loader.Archs4ExpressionDataLoaderFactory;
import org.reactome.idg.loader.GenePairCorrelationCalculator;

public class CorrelationVerification
{
	private static String pathToFile;
	private static final Logger logger = LogManager.getLogger();
	
	public static void main(String[] args) throws IOException
	{
		pathToFile = args[0];
		int numPairs = Integer.parseInt(args[1]);
		// First, we need to read the HDF file.
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(pathToFile);
		loader.loadMetaData();
		
		logger.info("Number of correlations to calculate: {}", numPairs);
		
		Random random = new Random(); // should I bother seeding this with the current time?
		// Need to know how large the random numbers can be.
		int maxNum = loader.getGeneIndices().keySet().size();
		// get a bunch of random numbers
		for (int i : random.ints(numPairs, 0, maxNum).toArray())
		{
			String geneSymbol1 = loader.getGeneIndicesToNames().get(i);
			String geneSymbol2 = geneSymbol1;
			int otherIndex = i;
			// Make sure we don't select the same two genes, or the correlation will just be 1.0, and that's not very interesting.
			while (geneSymbol2 == geneSymbol1)
			{
				otherIndex = random.nextInt(maxNum);
				geneSymbol2 = loader.getGeneIndicesToNames().get(otherIndex);
			}
			logger.info("Two genes selected: {} and {}", geneSymbol1, geneSymbol2);
			// Now that we have two genes, compute cross-tissue correlations.
			GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator(geneSymbol1, geneSymbol2, loader);
			double correlation = calculator.calculateGenePairCorrelation();
			logger.info("Correlation between {} and {} is: {}", geneSymbol1, geneSymbol2, correlation);
		}
	}
}
