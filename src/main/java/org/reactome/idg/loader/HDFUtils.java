package org.reactome.idg.loader;

import hdf.hdf5lib.H5;
import hdf.hdf5lib.HDF5Constants;

public class HDFUtils
{

	private HDFUtils() {}
	
	/**
	 * Reads data from a dataset.
	 * @param dataset_id - The dataset ID.
	 * @param space_id - The ID of the dataset's dataspace.
	 * @param dimx - The size of the x-dimension of the portion to read.
	 * @param dimy - The size of the y-dimension of the portion to read.
	 * @return An array that is <code>dimx</code> x <code>dimy</code>. 
	 */
	public static double[][] readData(long dataset_id, long space_id, int dimx, int dimy)
	{
		// TODO: Maybe make an overload of this function where one array dimension is always 1 so that we can return a 1-dimensional array?
		double[][] dset_data = new double[dimx][dimy];
		int rank = 2;
		long[] dims = new long[rank];
		dims[0] = dimx;
		dims[1] = dimy;
		long memspace_id = H5.H5Screate_simple(rank, dims, null);
		long type_id = H5.H5Dget_type(dataset_id);
		H5.H5Dread(dataset_id, type_id, memspace_id, space_id, HDF5Constants.H5P_DEFAULT, dset_data);
		H5.H5Sclose(memspace_id);
		H5.H5Tclose(type_id);
		return dset_data;
	}
	
	/**
	 * Reads data from a dataset.
	 * @param dataset_id - The dataset ID.
	 * @param space_id - The ID of the dataset's dataspace.
	 * @param dimx - The size of the x-dimension of the portion to read.
	 * @return An array that has length of <code>dimx</code>. 
	 */
	public static double[] readData(long dataset_id, long space_id, int dimx)
	{
		double[] dset_data = new double[dimx];
		int rank = 2;
		long[] dims = new long[rank];
		dims[0] = dimx;
		dims[1] = 1;
		long memspace_id = H5.H5Screate_simple(rank, dims, null);
		long type_id = H5.H5Dget_type(dataset_id);
		H5.H5Dread(dataset_id, type_id, memspace_id, space_id, HDF5Constants.H5P_DEFAULT, dset_data);
		H5.H5Sclose(memspace_id);
		H5.H5Tclose(type_id);
		return dset_data;
	}
	
	/**
	 * Reads an array from dataset into an array of StringBuffers.
	 * @param dsName - The name of the dataset in the HDF5 file. For example: "/meta/genes".
	 * @param dsSize - The number of data elements to read. This should be set to the size of the dataset.
	 * @return An array of StringBuffers.
	 */
	public static StringBuffer[] readDataSet(String fileName, String dsName, int dsSize)
	{
		long file_id = H5.H5Fopen(fileName, HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
		long dataset_id = H5.H5Dopen(file_id, dsName, HDF5Constants.H5P_DEFAULT);
		long type_id = H5.H5Dget_type(dataset_id);
		long dataWidth = H5.H5Tget_size(type_id);
		long space_id = H5.H5Dget_space(dataset_id);
		long[] dims = new long[2];
		long[] maxdims = new long[2];
		H5.H5Sget_simple_extent_dims(space_id, dims, maxdims);
		byte[][] dset_data = new byte[dsSize][(int) dataWidth];
		StringBuffer[] str_data = new StringBuffer[(int) maxdims[0]];
		H5.H5Dread(dataset_id, type_id, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, dset_data);
		byte[] tempbuf = new byte[(int) dataWidth];
		for (int indx = 0; indx < dset_data.length; indx++)
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
		return str_data;
	}
}
