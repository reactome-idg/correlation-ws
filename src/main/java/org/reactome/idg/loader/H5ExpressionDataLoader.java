package org.reactome.idg.loader;

import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;

public class H5ExpressionDataLoader
{
	private static final Logger logger = LogManager.getLogger();
	// Need a list of all Tissue-type names.
	private Set<String> tissueTypes = new HashSet<>();
	// A mapping of gene name to row-index in expression dataset.
	private Map<String, Integer> geneIndices = new HashMap<>();
	// A mapping of indices and the tissue type associated with it.
	private Map<Integer, String> indexOfTissues = new HashMap<>();
	private Map<String, List<Integer>> tissueTypeToIndex = new HashMap<>();

	// TODO: NUM_ROWS can probably be queried directly from the file, no need to hard-code it.
	private static final int NUM_ROWS=167726;
	private static String FILENAME = "/media/sshorser/data/reactome/IDGFiles/human_matrix.h5";
	String[] Xname = new String[1];
	// Some data-set names we will be using.
	private final static String expressionDSName = "/data/expression";
	private final static String genesDSName = "/meta/genes";
	private final static String tissueDSName = "/meta/Sample_source_name_ch1";
	// TODO: NUM_GENES_IN_FILE can probably be queried from the file. No need to hard code it.
	private static final int NUM_GENES_IN_FILE = 35238;
	
	public int[][] getExpressionValuesforTissue(String tissue)
	{
		List<Integer> indicesForTissue = this.tissueTypeToIndex.get(tissue);
		logger.info("number of samples for tissue ({}): {}", tissue, indicesForTissue.size());
		int[][] expressionValues = new int[NUM_GENES_IN_FILE][indicesForTissue.size()];
		
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
		long dset_space_id = H5.H5Dget_space(dataset_id);
		
		// Extract the elements for each gene. Since the gaps between the tissue columns is probably not regular enough to work with hyperslabs, we'll just do all
		// elements for each gene.
		int geneCount = 0;
		for (String gene : geneIndices.keySet())
		{
//			List<Integer> expressionValuesForGene = getExpressionValuesForGeneAndTissue(gene, tissue);
			int[] expressionValuesForGene = getExpressionValuesForGeneAndTissue(gene, tissue);
			geneCount++;
//			for (int i = 0; i < expressionValuesForGene.length; i++)
//			{
//				expressionValues[geneIndices.get(gene)][i] = expressionValuesForGene[i];
//			}
			if (geneCount % 100 == 0)
			{
				logger.debug("{} genes processed.", geneCount);
			}
				
			expressionValues[geneIndices.get(gene)] = expressionValuesForGene;
		}
		return expressionValues;
	}
	
	/**
	 * Gets a list of expression values for a gene/tissue pair. 
	 * @param gene
	 * @param tissue
	 * @return
	 */
	public int[] getExpressionValuesForGeneAndTissue(String gene, String tissue)
	{
		
		
		int geneIndex = geneIndices.get(gene);
		
//		String tissue = indexOfTissues.get(geneIndex);

		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
//		long type_id = H5.H5Dget_type(dataset_id);
//		long dataWidth = H5.H5Tget_size(type_id);
		
		
		List<Integer> indicesForTissue = this.tissueTypeToIndex.get(tissue);
//		int xoffset = indicesForTissue.get(0);
		
//		long[] start = { xoffset, geneIndex };
//		long[] stride = { 1, 1 };
//		long[] count = { 1, 1 };
//		long[] block = { indicesForTissue.get(indicesForTissue.size()-1) , 1};
		
		// this will create a hyperslab from the first occurrence of the tissue type to the last. BUT there might be some in the middle that are NOT
		// for this tissue, so we will need to filter it further, after the data is loaded.
		long space_id = H5.H5Dget_space(dataset_id);
//		int status = H5.H5Sselect_hyperslab(space_id, HDF5Constants.H5S_SELECT_SET, start, null, count, block);
		long[][] elementCoords = new long[indicesForTissue.size()][2];
		
		for (int i = 0; i < indicesForTissue.size(); i++)
		{
			elementCoords[i][1] = geneIndex;
			elementCoords[i][0] = indicesForTissue.get(i);
		}
		
		int status = H5.H5Sselect_elements(space_id, HDF5Constants.H5S_SELECT_SET, indicesForTissue.size(), elementCoords);
//		System.out.println("Is selection valid? " + H5.H5Sselect_valid(space_id));
		if (status < 0)
		{
			logger.error("Error selecting elements! response code: {}",status);
			return null;
		}
		int dimx = indicesForTissue.size();
		int dimy = 1;
		long type_id = H5.H5Dget_type(dataset_id);
//		long dataWidth = H5.H5Tget_size(type_id);
		int[][] dset_data = new int[dimx][dimy];
		
//		long memtype_id = H5.H5Tcopy(HDF5Constants.H5T_NATIVE_INT);
//		long memtype_id = H5.H5Tcopy(type_id);
//		
//		H5.H5Tset_size(memtype_id, dataWidth);
		long[] dims = new long[2];
		dims[0] = dimx;
		dims[1] = 1;
		long memspace_id = H5.H5Screate_simple(2, dims, dims);
		H5.H5Dread(dataset_id, type_id, memspace_id, space_id, HDF5Constants.H5P_DEFAULT, dset_data);
		int[] expressionValues = new int[dimx];
		for (int i = 0; i < dimx; i ++)
		{
			for (int j = 0; j < dimy; j++)
			{
				int expressionValue = dset_data[i][j];
				// If the index is in the tissue type's list, then the expression value can go into the output.
//				if (indicesForTissue.contains(expressionValue))
				{
					expressionValues[i] = expressionValue;
//					expressionValues.add(expressionValue);
				}
			}
		}
		H5.H5close();
		return expressionValues;
	}
	
	public void loadGeneNames()
	{
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, genesDSName, HDF5Constants.H5P_DEFAULT);
		long type_id = H5.H5Dget_type(dataset_id);
		long dataWidth = H5.H5Tget_size(type_id);
		long space_id = H5.H5Dget_space(dataset_id);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		byte[][] dset_data = new byte[NUM_ROWS][(int) dataWidth];
		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
		byte[] tempbuf = new byte[(int) dataWidth];
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			for (int jndx = 0; jndx < dataWidth; jndx++)
			{
				tempbuf[jndx] = dset_data[indx][jndx];
			}
			str_data[indx] = new StringBuffer(new String(tempbuf).trim());
		}
		H5.H5Dclose(dataset_id);
		H5.H5Tclose(type_id);
		H5.H5Sclose(space_id);
		H5.H5Fclose(file_id);
		logger.info("Number of elements: {}",maxdims[0]);
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			this.geneIndices.put(str_data[indx].toString(), indx);
		}
		logger.info("Number of genes loaded: {}", this.geneIndices.size());
	}
	
	public void loadTissueTypeNames()
	{
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, tissueDSName, HDF5Constants.H5P_DEFAULT);
		long type_id = H5.H5Dget_type(dataset_id);
		long dataWidth = H5.H5Tget_size(type_id);
		long space_id = H5.H5Dget_space(dataset_id);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		byte[][] dset_data = new byte[NUM_ROWS][(int) dataWidth];
		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
		byte[] tempbuf = new byte[(int) dataWidth];
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			for (int jndx = 0; jndx < dataWidth; jndx++)
			{
				tempbuf[jndx] = dset_data[indx][jndx];
			}
			str_data[indx] = new StringBuffer(new String(tempbuf).trim());
		}
		H5.H5Dclose(dataset_id);
		H5.H5Tclose(type_id);
		H5.H5Sclose(space_id);
		H5.H5Fclose(file_id);
		logger.info("Number of elements: ", maxdims[0]);
		for (int indx = 0; indx < maxdims[0]; indx++)
		{
			String tissueType = str_data[indx].toString();
			this.tissueTypes.add(tissueType);
			this.indexOfTissues.put(indx, tissueType);
			if (this.tissueTypeToIndex.containsKey(tissueType))
			{
				this.tissueTypeToIndex.get(tissueType).add(indx);
			}
			else
			{
				List<Integer> list = new ArrayList<>();
				list.add(indx);
				this.tissueTypeToIndex.put(tissueType, list);
			}
		}
		logger.info("Number of distinct elements: {}", this.tissueTypes.size());
	}
	
	public void getExpressionValuesForGene(String geneId)
	{
		long file_id = H5.H5Fopen(FILENAME, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		
		if (file_id >= 0)
		{
			long dataset_id = H5.H5Dopen(file_id, genesDSName, HDF5Constants.H5P_DEFAULT);
			
			long type_id = H5.H5Dget_type(dataset_id);
			System.out.println(type_id);
			long dataWidth = H5.H5Tget_size(type_id);
			long space_id = H5.H5Dget_space(dataset_id);
			long[] dims = new long[2];
			long[] maxdims = new long[2];
			H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
			System.out.println(dims[0]);
			System.out.println(maxdims[0]);
			byte[][] dset_data = new byte[NUM_ROWS][(int) dataWidth];
			System.out.println(dataWidth);
			if (dataset_id >= 0)
			{
				long dcpl_id = H5.H5Dget_create_plist(dataset_id);
				if (dcpl_id >= 0)
				{
					if (dataset_id >= 0)
					{
						StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
						long memtype_id = H5.H5Tcopy(HDF5Constants.H5T_C_S1);
						H5.H5Tset_size(memtype_id, dataWidth);
						H5.H5Dread(dataset_id, memtype_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
						
						byte[] tempbuf = new byte[(int) dataWidth];
						for (int indx = 0; indx < maxdims[0]; indx++)
						{
							for (int jndx = 0; jndx < dataWidth; jndx++)
							{
								tempbuf[jndx] = dset_data[indx][jndx];
							}
							str_data[indx] = new StringBuffer(new String(tempbuf).trim());
						}
						
						for (int indx = 0; indx < maxdims[0]; indx++)
						{
							logger.debug("{} [{}]: {}", genesDSName, indx, str_data[indx]);
						}
//						System.out.println();
					}
				}
			}
			
			dataset_id = H5.H5Dopen(file_id, expressionDSName, HDF5Constants.H5P_DEFAULT);
			type_id = H5.H5Dget_type(dataset_id);
			System.out.println(type_id);
			System.out.println(H5.H5Tget_tag(type_id));
		}
	}

}
