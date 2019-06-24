package org.reactome.idg;

import java.io.IOException;

import org.junit.Test;
import org.reactome.idg.loader.GenePairCorrelationMatrixCalculator;

@SuppressWarnings("static-method")
public class TestCorrelationCalculator
{

	@Test
	public void testCorrelationCalculationsIT() throws IOException
	{
		// It takes about 20 mintues to calculate all gene pair-wise correlations for the heart.txt samples, 141 samples for all genes (~35k) in the H5 file.
		GenePairCorrelationMatrixCalculator heartCorMatrixCalculator = new GenePairCorrelationMatrixCalculator("src/test/resources/heart.txt", "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5");
		
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
	
}
