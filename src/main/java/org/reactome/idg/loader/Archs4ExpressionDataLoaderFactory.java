package org.reactome.idg.loader;

import java.util.HashMap;
import java.util.Map;

/**
 * Creates instances of Archs4ExpressionDataLoader
 * @author sshorser
 *
 */
public class Archs4ExpressionDataLoaderFactory
{
	// there could be multiple files, but there should only be ONE loader per file.
	private static Map<String, Archs4ExpressionDataLoader> loaders = new HashMap<>();
	
	private Archs4ExpressionDataLoaderFactory()
	{
		// private constructor.
	}
	
	/**
	 * Creates a new data loader for an HDF file. If a loader for this file
	 * has already been created, then that loader will be returned. If no loader
	 * has been created, then a new loader will be created and returned.
	 * @param pathToFile
	 * @return
	 */
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
