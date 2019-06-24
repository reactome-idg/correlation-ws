package org.reactome.idg;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;
import org.reactome.idg.loader.Archs4ExpressionDataLoader;
import org.reactome.idg.loader.Archs4ExpressionDataLoaderFactory;
import org.reactome.idg.loader.GenePairCorrelationCalculator;
import org.reactome.idg.loader.GenePairCorrelationMatrixCalculator;

@SuppressWarnings("static-method")
public class TestCorrelationCalculator
{

	private static final String PATH_TO_HDF = "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5";
	@Test
	public void testCorrelationCalculationsIT() throws IOException
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		// It takes about 20 mintues to calculate all gene pair-wise correlations for the heart.txt samples, 141 samples for all genes (~35k) in the H5 file.
		GenePairCorrelationMatrixCalculator heartCorMatrixCalculator = new GenePairCorrelationMatrixCalculator("src/test/resources/heart.txt", loader);
		
		double[][] d = heartCorMatrixCalculator.calculateCorrelation();
		for (int i = 0; i < Math.min(d.length, 20); i++)
		{
			for (int j = 0; j < Math.min(d[i].length, 20); j++)
			{
				System.out.print(d[i][j] + "\t");
			}
			System.out.println();
		}
	}
	
	@Test
	public void testGenePairCorrelationsIT() throws IOException
	{
		Archs4ExpressionDataLoader loader = Archs4ExpressionDataLoaderFactory.buildInstanceForHDFFile(PATH_TO_HDF);
		GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator("src/test/resources/heart.txt", "A1BG","A1BG", loader);
		double d = calculator.calculateGenePairCorrelation();
		System.out.println(d);
		
		assertEquals(d, 1.0, 0);
		
		d = calculator.calculateGenePairCorrelation("A1BG","A1CF");
		System.out.println(d);
		
		try
		{
			d = calculator.calculateGenePairCorrelation("BLAHBLAHBLAH","A1CF");
			System.out.println(d);
		}
		catch (IllegalArgumentException e)
		{
			e.printStackTrace();
			assertTrue(e.getMessage().contains("Gene BLAHBLAHBLAH is not recognized in the HDF file."));
		}
		
		try
		{
			d = calculator.calculateGenePairCorrelation("A1CF", "BLAH6666666");
			System.out.println(d);
		}
		catch (IllegalArgumentException e)
		{
			e.printStackTrace();
			assertTrue(e.getMessage().contains("Gene BLAH6666666 is not recognized in the HDF file."));
		}
	}
}
