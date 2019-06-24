package org.reactome.idg.loader;

import java.util.HashMap;
import java.util.Map;

public class Archs4ExpressionDataLoaderFactory
{
	// there could be multiple files, but there should only be ONE loader per file.
	private static Map<String, Archs4ExpressionDataLoader> loaders = new HashMap<>();
	
	private Archs4ExpressionDataLoaderFactory()
	{
		// private constructor.
	}
	
	public static Archs4ExpressionDataLoader buildInstanceForHDFFile(String pathToFile)
	{
		if (loaders.containsKey(pathToFile))
		{
			return loaders.get(pathToFile);
		}
		else
		{
			Archs4ExpressionDataLoader loader = new Archs4ExpressionDataLoader(pathToFile);
			loader.loadCounts();
			loader.loadMetaData();
			loaders.put(pathToFile, loader);
			return loader;
		}
	}
}
