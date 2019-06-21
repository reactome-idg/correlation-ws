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

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.junit.BeforeClass;
import org.junit.Test;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;

@SuppressWarnings("static-method")
public class TestH5Loader
{
	@BeforeClass
	public static void setup()
	{
		Archs4ExpressionDataLoader.setHdfExpressionFileAndLoadCounts("/media/sshorser/data/reactome/IDGFiles/human_matrix.h5");
	}
	
//	@Test
//	public void testGetExpressionsForGene()
//	{
//		Archs4ExpressionDataLoader.getExpressionValuesForGene("AAA");
//	}
	
	@Test
	public void testH5LoaderTissueTypesIT()
	{
		Archs4ExpressionDataLoader.loadMetaData();
		assertNotNull(Archs4ExpressionDataLoader.getTissueTypes());
		assertNotNull(Archs4ExpressionDataLoader.getTissueTypeToIndex());
	}
	
	@Test
	public void testH5LoaderGeneNamesIT()
	{
		Archs4ExpressionDataLoader.loadMetaData();
		assertNotNull(Archs4ExpressionDataLoader.getGeneIndices());
	}
	
	@Test
	public void testH5LoaderGetExpressionsForGeneAndTissue()
	{
		Archs4ExpressionDataLoader.loadMetaData();
		int expressionValues[] = Archs4ExpressionDataLoader.getExpressionValuesForGeneAndTissue("A1BG", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(124, expressionValues[0]);
		System.out.println("\n");
		
		
		expressionValues = Archs4ExpressionDataLoader.getExpressionValuesForGeneAndTissue("A2MP1", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(119,  expressionValues[0]);
		System.out.println("\n");

		expressionValues = Archs4ExpressionDataLoader.getExpressionValuesForGeneAndTissue("A1BG", "Heart");
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
		Archs4ExpressionDataLoader.loadMetaData();

		String tissueName = "HeLa ELAVL1/HuR siRNA1 5d";
		loadTissueType(tissueName);
		
		tissueName = "HeLa mock knockdown 5d";
		loadTissueType(tissueName);
		
		tissueName = "HSTL_Spleen";
		loadTissueType(tissueName);
		
		tissueName = "brain";
		loadTissueType(tissueName);
		
//		tissueName = "LCL-derived iPSC"; //biggest data set! takes > 1 hour on my workstation.
//		loadTissueType(tissueName);
	}

	private void loadTissueType(String tissueName)
	{
		LocalDateTime start;
		LocalDateTime end;
		int[][] expressionValues;
		start = LocalDateTime.now();
		expressionValues = Archs4ExpressionDataLoader.getExpressionValuesforTissue(tissueName);
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
	public void testCalculateCorrelationIT()
	{
		Archs4ExpressionDataLoader.loadMetaData();
		LocalDateTime start;
		LocalDateTime end;
		int[][] expressionValues;
		start = LocalDateTime.now();
		String tissueName = "HSTL_Spleen";
		expressionValues = Archs4ExpressionDataLoader.getExpressionValuesforTissue(tissueName);
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
		
		double[][] dExpressionValues = new double[expressionValues.length][expressionValues[0].length];
		for (int i = 0; i < expressionValues.length; i++)
		{
			for (int j = 0; j < expressionValues[i].length; j++)
			{
				dExpressionValues[i][j] = expressionValues[i][j];
			}
		}
		
		// this really just calculates the correlation for a single tissue, and I'm not 100% sure this is the correct way to do it.
		start = LocalDateTime.now();
		PearsonsCorrelation cor = new PearsonsCorrelation(dExpressionValues);
		RealMatrix corMatrix = cor.getCorrelationMatrix();
		end = LocalDateTime.now();
		System.out.println("Elapsed time: " + Duration.between(start, end).toString());
		System.out.println("Column dimensions: " + corMatrix.getColumnDimension());
		System.out.println("Row dimensions: " + corMatrix.getRowDimension());
		
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				System.out.print(corMatrix.getEntry(i, j)+"\t\t");
			}
			System.out.print("\n");
		}
	}
	
	@Test
	public void testTissueTypes()
	{
		Archs4ExpressionDataLoader.loadMetaData();
		Archs4ExpressionDataLoader.getTissueTypes().stream().sorted().collect(Collectors.toList()).forEach(t -> System.out.println(t));
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
		Archs4ExpressionDataLoader.loadMetaData();
		Map<String, List<Integer>> tissueMap = Archs4ExpressionDataLoader.getTissueTypeToIndex();
		tissueMap.keySet().stream().sorted().forEach(t -> System.out.println(t + "\t" + tissueMap.get(t).size()));
	}
	
	@Test
	public void getSampleIDsIT()
	{
		Archs4ExpressionDataLoader.loadSampleIndices();
	}
	
	@Test
	public void getExpressionValuesForTissueFromSampleIDFile() throws IOException
	{
		Archs4ExpressionDataLoader.loadMetaData();
		int[][] expressionValues = Archs4ExpressionDataLoader.getExpressionValuesforTissue(Paths.get("src/test/resources/heart.txt"));
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
}
