package org.reactome.idg;

import java.io.IOException;

import org.junit.Test;
import org.reactome.idg.loader.GenePairCorrelationCalculator;

@SuppressWarnings("static-method")
public class TestCorrelationCalculator
{

	@Test
	public void testCorrelationCalculationsIT() throws IOException
	{
		GenePairCorrelationCalculator calculator = new GenePairCorrelationCalculator("A1BG", "A1CF", "src/test/resources/heart.txt", "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5");
		
		double[][] d = calculator.calculateCorrelation();
		for (int i = 0; i < Math.min(d.length, 20); i++)
		{
			for (int j = 0; j < Math.min(d[i].length, 20); i++)
			{
				System.out.print(d[i][j] + "\t");
			}
			System.out.println();
		}
	}
	
}
