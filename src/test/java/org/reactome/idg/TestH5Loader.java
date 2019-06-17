package org.reactome.idg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.reactome.idg.loader.H5ExpressionDataLoader;

@SuppressWarnings("static-method")
public class TestH5Loader
{
	@Test
	public void testH5LoaderIT()
	{
		H5ExpressionDataLoader loader = new H5ExpressionDataLoader();
		loader.getExpressionValuesForGene("AAA");
	}
	
	@Test
	public void testH5LoaderTissueTypesIT()
	{
		H5ExpressionDataLoader loader = new H5ExpressionDataLoader();
		loader.loadTissueTypeNames();
	}
	
	@Test
	public void testH5LoaderGeneNamesIT()
	{
		H5ExpressionDataLoader loader = new H5ExpressionDataLoader();
		loader.loadGeneNames();
	}
	
	@Test
	public void testH5LoaderGetExpressionsForGeneAndTissue()
	{
		H5ExpressionDataLoader loader = new H5ExpressionDataLoader();
		loader.loadGeneNames();
		loader.loadTissueTypeNames();
		int expressionValues[] = loader.getExpressionValuesForGeneAndTissue("A1BG", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(124, expressionValues[0]);
		System.out.println("\n");
		
		
		expressionValues = loader.getExpressionValuesForGeneAndTissue("A2MP1", "HeLa ELAVL1/HuR siRNA1 5d");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.length > 0);
		System.out.println(expressionValues.length);
		Arrays.stream(expressionValues).forEach(x -> System.out.print(x+", "));
		assertEquals(119,  expressionValues[0]);
		System.out.println("\n");
/*
		expressionValues = loader.getExpressionValuesForGeneAndTissue("A1BG", "Heart");
		assertNotNull(expressionValues);
		assertTrue(expressionValues.size() > 0);
		System.out.println(expressionValues.size());
		expressionValues.stream().forEach(x -> System.out.print(x+", "));
		System.out.println("\n");
		
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
		H5ExpressionDataLoader loader = new H5ExpressionDataLoader();
		H5ExpressionDataLoader.loadGeneNames();
		H5ExpressionDataLoader.loadTissueTypeNames();
//		LocalDateTime start, end;
//		int[][] expressionValues;
//		String tissueName = "LCL-derived iPSC"; //biggest data set!
//		String tissueName = "HSTL_Spleen";
		String tissueName = "HeLa ELAVL1/HuR siRNA1 5d";
		loadTissueType(loader, tissueName);
		
		tissueName = "HeLa mock knockdown 5d";
		loadTissueType(loader, tissueName);
		
		tissueName = "HSTL_Spleen";
		loadTissueType(loader, tissueName);
		
		tissueName = "brain";
		loadTissueType(loader, tissueName);
		
		tissueName = "LCL-derived iPSC";
		loadTissueType(loader, tissueName);
	}

	private void loadTissueType(H5ExpressionDataLoader loader, String tissueName)
	{
		LocalDateTime start;
		LocalDateTime end;
		int[][] expressionValues;
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
}
