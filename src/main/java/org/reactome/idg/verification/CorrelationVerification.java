package org.reactome.idg.verification;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;

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
		Map<String, Double> calculatedCorrelations = new HashMap<>(numPairs);
		// First, we need to read the HDF file.
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(pathToHDFFile);
		loader.loadMetaData();
		
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
		for (int i : random.ints(numPairs, 0, maxNum).toArray())
		{
			int geneIndex1 = 0;//i;
			String geneSymbol1 = loader.getGeneIndicesToNames().get(geneIndex1);
			while (!geneNamesToIndices.containsKey(geneSymbol1))
			{
				geneIndex1 = random.nextInt(maxNum);
				geneSymbol1 = loader.getGeneIndicesToNames().get(geneIndex1);
			}
			String geneSymbol2 = "A1CF";//geneSymbol1;
			int otherIndex = 1;//geneIndex1;
			// Make sure we don't select the same two genes, or the correlation will just be 1.0, and that's not very interesting.
//			while (geneSymbol2 == geneSymbol1 || !geneNamesToIndices.containsKey(geneSymbol2))
//			{
//				otherIndex = random.nextInt(maxNum);
//				geneSymbol2 = loader.getGeneIndicesToNames().get(otherIndex);
//			}
			logger.info("Two genes selected: {} and {}", geneSymbol1, geneSymbol2);
			// Now that we have two genes, compute cross-tissue correlations.
			GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator(geneSymbol1, geneSymbol2, loader);
			double correlation = calculator.calculateGenePairCorrelation();
			logger.info("Correlation between {} and {} is: {}", geneSymbol1, geneSymbol2, correlation);
			calculatedCorrelations.put(DataRepository.generateKey(geneSymbol1, geneSymbol2), correlation);
		}
		
		
		Map<String,Double> correlationsFromFile = new HashMap<>(numPairs);
		for (String geneKey : calculatedCorrelations.keySet())
		{
			String[] parts = geneKey.split("\\|");
			String gene1 = parts[0];
			String gene2 = parts[1];
			int geneIndex1 = geneNamesToIndices.get(gene1);
			int geneIndex2 = geneNamesToIndices.get(gene2);
			int lineCount = 0;
			// now that we have the indices, we need to read into the file.
			// need to find a way to read to a specific line number in a file, using Java.
			try(FileInputStream fis = new FileInputStream(pathToCorrelationFile);
					Scanner scanner = new Scanner(fis);)
			{
				boolean done = false;
				while (scanner.hasNext() && !done)
				{
					String line = scanner.nextLine();
					if (line.startsWith("\""+gene1+"\""))
					{
						// now we found the line we need, get the correlation value.
						String[] lineParts = line.split(",");
						double correlationFromFile = Double.parseDouble(lineParts[geneIndex2]);
						logger.info("Correlation (according to the FILE) between {} and {} is: {}", gene1, gene2, correlationFromFile);
						done = true;
						correlationsFromFile.put(DataRepository.generateKey(gene1, gene2), correlationFromFile);
						logger.info("Delta: {}", calculatedCorrelations.get(DataRepository.generateKey(gene1, gene2)) - correlationsFromFile.get(DataRepository.generateKey(gene1, gene2)) );
					}
//					lineCount++;
				}
				
			}
		}
	}
}
