package org.reactome.idg.model;

import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Table;

/**
 * Represents the provenance of a gene-pair correlation.
 * Attributes of Provenance are name, URL, category, subcategory.
 * @author sshorser
 *
 */
@Entity
@Table(name = "provenance"/*
							 * , indexes = { @Index(columnList = "name"), @Index(columnList =
							 * "url"), @Index(columnList = "category") }
							 */)
public class Provenance implements Serializable
{
	/**
	 * generated serialVersionUID
	 */
	private static final long serialVersionUID = 2045677264758847791L;

	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private long id;
	
	@Column(nullable = false)
	private String name;
	@Column
	private String url;
	
	// Should these be free-form text fields or should they be foreign keys to proper classes?
	@Column
	private String category;
	@Column
	private String subcategory;
	
//	@OneToMany(mappedBy = "provenance_id")
//	Set<GenePairCorrelation> correlations;
	
	public Provenance(String name, String URL, String category, String subcategory)
	{
		this.setName(name);
		this.setUrl(URL);
		this.setCategory(category);
		this.setSubcategory(subcategory);
	}

	public Provenance()
	{
	}

	public String getName()
	{
		return name;
	}
	
	public void setName(String name)
	{
		this.name = name;
	}
	public String getUrl()
	{
		return url;
	}
	
	public void setUrl(String url)
	{
		this.url = url;
	}
	
	public String getCategory()
	{
		return category;
	}
	
	public void setCategory(String category)
	{
		this.category = category;
	}
	
	public String getSubcategory()
	{
		return subcategory;
	}
	
	public void setSubcategory(String subcategory)
	{
		this.subcategory = subcategory;
	}
	
	public long getId()
	{
		return id;
	}

	public void setId(long id)
	{
		this.id = id;
	}
	
	//TODO: Add: species, tissue-type, other fields...
}
