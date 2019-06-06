package org.reactome.idg;

import java.util.Optional;

/**
 * Represents the provenance of a dataset.
 * Attributeas include: Name of dataset, URL where found, Category, Subcategory.
 * @author sshorser
 * @deprecated Use {@link org.reactome.idg.model.Provenance}
 * @param <T>
 *
 */
public class Provenance implements Comparable<Provenance>
{
	private static final String TOKEN_FOR_NULL_VALUES = "N/A";
	private String name;
	private String URL;
	private String category;
	private String subcategory;
	// TODO: Add other fields for things such as species and tissue type.
	
	/**
	 * Creates a new Provenance object. NULL values will be replaced with "N/A"
	 * @param name
	 * @param URL
	 * @param category
	 * @param subcategory
	 */
	public Provenance(String name, String URL, String category, String subcategory)
	{
		this.name = Optional.ofNullable(name).orElse(TOKEN_FOR_NULL_VALUES);
		this.URL = Optional.ofNullable(URL).orElse(TOKEN_FOR_NULL_VALUES);
		this.category = Optional.ofNullable(category).orElse(TOKEN_FOR_NULL_VALUES);
		this.subcategory = Optional.ofNullable(subcategory).orElse(TOKEN_FOR_NULL_VALUES);
	}
	
	public String getName()
	{
		return name;
	}

	public String getURL()
	{
		return URL;
	}

	public String getCategory()
	{
		return category;
	}

	public String getSubcategory()
	{
		return subcategory;
	}

	/**
	 * The String representation of a Provenance object is: "Name: $name ; URL: $url ; Category: $category ; Subcategory: $subcategory ;" 
	 */
	@Override
	public String toString()
	{
		return "Name: "+this.name + " ; URL: " + this.URL + " ; Category: " + this.category + " ; Subcategory: " + this.subcategory + " ;";
	}
	
	/**
	 * Compares two Provenance objects. They are compared by their "toString()" values.
	 * @param other - the Other provenance object.
	 * @return a negative value if "this" is > than other. 0 if they are the same. A positive value if "this" is < than the other.
	 */
	@Override
	public int compareTo(Provenance other)
	{

		int comparisonValue;

		comparisonValue = this.toString().compareTo(other.toString());
		
		return comparisonValue;
	}
	/*	
	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof Provenance)
		{
			return this.name.equals( ((Provenance)obj).getName() )
				&& this.URL.equals( ((Provenance)obj).getURL() )
				&& this.category.equals( ((Provenance)obj).getCategory() )
				&& this.subcategory.equals( ((Provenance)obj).getSubcategory());
		}
		return false;
	}
	

	@Override
	public int hashCode()
	{
		return Objects.hash(this.name, this.URL, this.category, this.subcategory);
	}
	*/
}
