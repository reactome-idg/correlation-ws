package org.reactome.idg.model;

import javax.persistence.Column;
import javax.persistence.ConstraintMode;
import javax.persistence.Entity;
import javax.persistence.ForeignKey;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.Index;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "gene_pair_correlation", indexes = { @Index(columnList = "gene_1, gene_2", unique = false),
													@Index(columnList = "gene_1,gene_2,provenance_id", unique = true, name = "idx_gene_pair_provenance")})
public class GenePairCorrelation
{
	public GenePairCorrelation()
	{
		// Hibernate requires a deafult no-arg constructor.
	}
	
	public GenePairCorrelation(String gene1, String gene2, double correlationValue, Provenance provenance)
	{
		this.gene1 = gene1;
		this.gene2 = gene2;
		this.correlationValue = correlationValue;
		this.provenance = provenance;
	}
	
	@Id
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id;

	@Column(name = "gene_1", nullable = false, length = 40)
	private String gene1;
	
	@Column(name = "gene_2", nullable = false, length = 40)
	private String gene2;
	
	@Column(name = "correlation_value", nullable = false, columnDefinition = "DOUBLE PRECISION(7,6)")
	private double correlationValue;
	
	@JoinColumn(name = "provenance_id", foreignKey = @ForeignKey(value = ConstraintMode.CONSTRAINT), nullable = false)
	@ManyToOne
	private Provenance provenance;

	public GenePairCorrelation(Long id)
	{
		this.id = id;
	}
 
	public Long getId()
	{
		return id;
	}

	public String getGene1()
	{
		return gene1;
	}

	public void setGene1(String gene1)
	{
		this.gene1 = gene1;
	}

	public String getGene2()
	{
		return gene2;
	}

	public void setGene2(String gene2)
	{
		this.gene2 = gene2;
	}

	public double getCorrelationValue()
	{
		return correlationValue;
	}

	public void setCorrelationValue(double correlationValue)
	{
		this.correlationValue = correlationValue;
	}

	public Provenance getProvenance()
	{
		return provenance;
	}

	public void setProvenance(Provenance provenance)
	{
		this.provenance = provenance;
	}
	
}