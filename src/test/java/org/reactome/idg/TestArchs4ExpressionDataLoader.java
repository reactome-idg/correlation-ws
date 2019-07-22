package org.reactome.idg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Paths;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.junit.BeforeClass;
import org.junit.Test;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;
import org.reactome.idg.loader.Archs4ExpressionDataLoaderFactory;

@SuppressWarnings("static-method")
public class TestArchs4ExpressionDataLoader
{
	final static String PATH_TO_HDF = "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5";
	
	@BeforeClass
	public static void setup()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
	}

	
	@Test
	public void testH5LoaderGetExpressionsForGeneAndTissue()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		double expressionValues[] = loader.getExpressionValuesForGeneAndTissue("A1BG", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(124, expressionValues[0], 0.00000000001);
		System.out.println("\n");
		
		
		expressionValues = loader.getExpressionValuesForGeneAndTissue("A2MP1", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(119,  expressionValues[0], 0.00000000001);
		System.out.println("\n");

		expressionValues = loader.getExpressionValuesForGeneAndTissue("A1BG", "Heart");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
/*		
		// Note: "heart" and "Heart" will yield different results. Tissue-type names are case-sensitive in the source file.
		expressionValues = loader.getExpressionValuesForGeneAndTissue("A1CF", "heart");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
		
		expressionValues = loader.getExpressionValuesForGeneAndTissue("A1CF", "brain");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
		
		expressionValues = loader.getExpressionValuesForGeneAndTissue("BRCA1", "colon");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
		
		expressionValues = loader.getExpressionValuesForGeneAndTissue("BRCA2", "colon");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");

		expressionValues = loader.getExpressionValuesForGeneAndTissue("A1BG", "colon");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
		*/
	}
	
	@Test
	public void testGetExpressionValuesForTissueAllGenes()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);

		String tissueName = "HeLa ELAVL1/HuR siRNA1 5d";
		loadTissueType(tissueName, loader);
		
		tissueName = "HeLa mock knockdown 5d";
		loadTissueType(tissueName, loader);
		
		tissueName = "HSTL_Spleen";
		loadTissueType(tissueName, loader);
		
		tissueName = "brain";
		loadTissueType(tissueName, loader);
		
//		tissueName = "LCL-derived iPSC"; //biggest data set! takes > 1 hour on my workstation.
//		loadTissueType(tissueName);
	}

	private void loadTissueType(String tissueName, Archs4ExpressionDataLoader loader)
	{
		LocalDateTime start;
		LocalDateTime end;
		double[][] expressionValues;
		start = LocalDateTime.now();
		expressionValues = loader.getExpressionValuesforTissue(tissueName);
		end = LocalDateTime.now();
		System.out.println("Elapsed time: " + Duration.between(start, end).toString());

		System.out.println("Size: " + expressionValues.length + " x " + expressionValues[0].length);
		for (int i = 0; i < Math.min(10,expressionValues.length); i++)
		{
			for (int j = 0; j < Math.min(10, expressionValues[i].length); j++)
			{
				System.out.print(expressionValues[i][j]+"\t");
			}
			System.out.print("\n");
		}
	}
	
	@Test
	public void testTissueTypes()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		loader.getTissueTypes().stream().sorted().collect(Collectors.toList()).forEach(t -> System.out.println(t));
	}
//	
//	@Test
//	public void testGetTissueSampleMapping()
//	{
//		Archs4ExpressionDataLoader.loadTissueTypeNames();
//		Archs4ExpressionDataLoader.loadGeneNames();
//		
//		Map<String, String> map = Archs4ExpressionDataLoader.getTissueToSamples();
//		for (String tissue : map.keySet().parallelStream().sorted().collect(Collectors.toList()))
//		{
//			System.out.println(tissue+"\t"+map.get(tissue));
//			
//		}
//	}
	
	@Test
	public void getTissuesWithCounts()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		Map<String, List<Integer>> tissueMap = loader.getTissueTypeToIndex();
		tissueMap.keySet().stream().sorted().forEach(t -> System.out.println(t + "\t" + tissueMap.get(t).size()));
	}
	
	@Test
	public void getSampleIDsIT()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		loader.loadSampleIndices();
	}
	
	@Test
	public void getExpressionValuesForTissueFromSampleIDFile() throws IOException
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		int[] sampleIndices = loader.getSampleIndicesForTissue(Paths.get("src/test/resources/heart.txt"));
		double[][] expressionValues = loader.getExpressionValuesforTissue(Paths.get("src/test/resources/heart.txt"));
		System.out.println("Size: " + expressionValues.length + " x " + expressionValues[0].length);

		System.out.print("\t\t");
		for (int i = 0; i < Math.min(expressionValues.length, 10); i++)
		{
			System.out.print(loader.getGeneIndicesToNames().get(i)+ "\t");
		}
		System.out.print("\n");
		
		
		for (int i = 0; i < Math.min(expressionValues.length, 10); i++)
		{
			System.out.print(loader.getSampleIndexToID().get(sampleIndices[i])+"\t");
			for (int j = 0; j < Math.min(10, expressionValues[i].length); j++)
			{
				System.out.print(expressionValues[i][j]+"\t");
			}
			System.out.print("\n");
		}
	}
	
	@Test
	public void testGetGeneIndicesToNamesIT()
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		assertTrue(loader.getGeneIndicesToNames().get(0).equals("A1BG"));
	}
}
