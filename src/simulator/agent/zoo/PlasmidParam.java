package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.agent.zoo.PlasmidBac;
import simulator.agent.Species;
import utils.LogFile;
import utils.UnitConverter;
import utils.XMLParser;

/**
 * \brief TODO
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidParam extends ActiveParam
{
	/**
	 * Info used in output writing.
	 */
	public String name;
	
	/**
	 * Default ID of this plasmid.
	 */
	public Integer plasmidID = 0;
	
	/**
	 * Default copy number of this plasmid in a host.
	 */
	public Integer copyNumDefault = 1;
	/**
	 * Type of cost of plasmid. Can be plasmidAdapt, hostAdapt, and constantCost.
	 * (Transferred from PlasmidBacParam to here, otherwise transconjugants do not have information about type of plasmid cost)
	 */
	public String costOfPlasmid;
	
	/** 
	 * The initial cost of the plasmid
	*/
	protected double initialCost;
	
	/**
	 * The minimal cost of a plasmid
	 */
	protected double basalCost;
	
	/**
	 * Rate at which the fitness cost of a plasmid decreases over time due to co-evolution with the host.
	 * When the plasmid enters the host for the first time, the fitness cost is equal to initial cost.
	 * Then fitness cost evolves, exponentially decreasing to the final, basal cost, with rate constant _rateCostDecrease
	 */
	protected double rateCostDecrease;
	
	/**
	 * Length of the pilus associated with this plasmid for conjugation 
	 * (in um). Note that the pilus should reach between cell surfaces, not
	 * between cell centres.
	 */
	public Double pilusLength = 5.0;
	
	/**
	 * After donating a plasmid of this species, a host needs to recover
	 * before it may exchange again. This parameter gives the time needed to
	 * recover (in hours).
	 * 
	 * TODO Note that this was previously called transferLag - protocol files
	 * will need to be modified.
	 */
	public Double donationLag = 5.0;
	
	/**
	 * After receiving a plasmid of this species, a host needs to recover
	 * before it may receive again. This parameter gives the time needed to
	 * recover (in hours).
	 */
	public Double receptionLag = 1.0;
	
	/**
	 * Probability that this plasmid will be lost by the host during division.
	 */
	public Double lossProbability = 0.0;
	
	/**
	 * Probability that a possible transfer will be successful.
	 */
	public Double transferProficiency = 1.0;
	
	/**
	 * Parameter for collision frequency for this plasmid in chemostat simulations
	 */
	public double collisionCoeff = 0.1;
	
/**
	 * The maximum scan speed (in cells/hour) for this plasmid in biofilm simulations
	 */
	public Double scanSpeed = 10.0;
	
	/**
	 * Plasmid markers (used in host compatibility).
	 */
	public ArrayList<String> hostCompatibilityMarkers =
													new ArrayList<String>();
	
	/**
	 * Plasmid Markers for plasmid compatibility (based on incompatibility
	 * groups).
	 */
	public ArrayList<String> plasmidCompatibilityMarkers =
													new ArrayList<String>();
	
	/**
	 * List of integer indices corresponding to the reactions that Plasmids
	 * of this species encode for.
	 */
	public ArrayList<Integer> reactionsEncoded = new ArrayList<Integer>();
	
	private static HashMap<Integer,String> plasmidIDAndNameList=new HashMap<>();
		
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidParam()
	{
		super();
	}
	
	/**
	 * Called during creation of the species
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		
		Integer tempInt;
		Double tempDbl;
		String tempStr;
		LinkedList<XMLParser> tempXML;
		
		name = aSpeciesRoot.getName();
		
		tempInt = getSpeciesParameterInteger("plasmidID",
				aSpeciesRoot, speciesDefaults);
		plasmidID = ( tempInt == null ) ? plasmidID : tempInt;
		plasmidIDAndNameList.put(plasmidID, name);		
		tempInt = getSpeciesParameterInteger("copyNumDefault",
											aSpeciesRoot, speciesDefaults);
		copyNumDefault = ( tempInt == null ) ? copyNumDefault : tempInt;
		
		tempDbl = getSpeciesParameterDouble("initialCost",
				aSpeciesRoot, speciesDefaults);
		tempStr = getSpeciesParameterString ("costOfPlasmid", aSpeciesRoot, speciesDefaults);
		tempStr = (tempStr==null) ? "Unknown" : tempStr;
		costOfPlasmid = tempStr;
		initialCost = Double.isFinite(tempDbl) ? tempDbl : initialCost;
		
		tempDbl = getSpeciesParameterDouble("basalCost",
				aSpeciesRoot, speciesDefaults);
		basalCost = Double.isFinite(tempDbl) ? tempDbl : basalCost;
				
		tempDbl = getSpeciesParameterDouble("rateCostDecrease",
				aSpeciesRoot, speciesDefaults);
		rateCostDecrease = Double.isFinite(tempDbl) ? tempDbl : rateCostDecrease;
		
		tempDbl = getSpeciesParameterLength("pilusLength", 
											aSpeciesRoot, speciesDefaults);
		pilusLength = Double.isFinite(tempDbl) ? tempDbl : pilusLength;
		
		tempDbl = getSpeciesParameterTime("donationLag",
											aSpeciesRoot, speciesDefaults);
		donationLag = Double.isFinite(tempDbl) ? tempDbl : donationLag;
		
		tempDbl = getSpeciesParameterTime("receptionLag",
											aSpeciesRoot, speciesDefaults);
		receptionLag = Double.isFinite(tempDbl) ? tempDbl : receptionLag;
		
		tempDbl = getSpeciesParameterDouble("lossProbability",
											aSpeciesRoot, speciesDefaults);
		lossProbability = Double.isFinite(tempDbl) ? tempDbl : lossProbability;

		tempDbl = getSpeciesParameterDouble("transferProficiency",
											aSpeciesRoot, speciesDefaults);
		transferProficiency = Double.isFinite(tempDbl) ?
												tempDbl : transferProficiency;
		tempDbl = getSpeciesParameterTime("scanSpeed",
											aSpeciesRoot, speciesDefaults);
		scanSpeed = Double.isFinite(tempDbl) ? tempDbl : scanSpeed;
		//collisionCoeff
		//LogFile.writeLog("collisionCoeff " + collisionCoeff);
		StringBuffer tempUnit = new StringBuffer("");
		tempDbl = aSpeciesRoot.getParamDbl("collisionCoeff", tempUnit);
		collisionCoeff = tempDbl*UnitConverter.time(tempUnit.toString());
		//LogFile.writeLog("collisionCoeff " + collisionCoeff);				
				
		/*
		 * 
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("reaction");
		for ( XMLParser parser : tempXML )
			if ( parser.getAttribute("status").equals("active") )
				reactionsEncoded.add(aSim.getReactionIndex(parser.getName()));
		/*
		 * If no host markers are given, assume all species may be hosts.
		 * Note: Old EpiBac protocol files should be unaffected by this.
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("compatibleHost");
		if ( tempXML == null || tempXML.isEmpty() )
		{
			//System.out.println("Entering all host compatibility markers");
			for ( Species spec : aSim.speciesList )
				if ( spec.getProgenitor() instanceof PlasmidBac )
					hostCompatibilityMarkers.add( spec.speciesName );
		}
		else
		{
			//System.out.println("Trying to read in host compatibility markers");
			for ( XMLParser parser : tempXML )
			{
				//System.out.println("\tAdding "+parser.getName());
				hostCompatibilityMarkers.add( parser.getName() );
			}
		}
		/*
		 * Try to add this plasmid name to the list of potential plasmids in 
		 * compatible hosts, if they have been initialised yet.
		 */
		for ( String hostName : hostCompatibilityMarkers )
			try
			{
				((PlasmidBacParam)aSim.getSpecies(hostName).getSpeciesParam())
											.potentialPlasmids.add(this.name);
			}
			catch ( Exception e) {}
		/*
		 * If no plasmid markers are given, assume all plasmids are
		 * incompatible.
		 * 
		 * TODO check that compatibility goes both ways?
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("compatiblePlasmid");
		if ( ! (tempXML == null || tempXML.isEmpty()) )
			for ( XMLParser parser : tempXML )
				plasmidCompatibilityMarkers.add( parser.getName() );
	}
	
	public boolean isHostCompatible(PlasmidBac targetRecipient)
	{
		//LogFile.writeLog("targetRecipient.getName() " + targetRecipient.getName());
		//LogFile.writeLog("hostCompatibilityMarkers " + hostCompatibilityMarkers.toString());
		if ( ! hostCompatibilityMarkers.contains(targetRecipient.getName()) )
			return false;
		return true;    
	}

	/**
	 * \brief Check if a potential host is compatible with a plasmid of this
	 * species.
	 * Update: Take Entry Exclusion System into consideration. If the same plasmid 
	 * already exists in the targetRecipient, then don't transfer. 
	 * 
	 * @param targetRecipient PlasmidBac to check compatibility with.
	 * @return boolean: true if target is compatible, false if not. 
	 */
	public boolean isCompatible(PlasmidBac targetRecipient)
	{
		if ( ! hostCompatibilityMarkers.contains(targetRecipient.getName()) )
			return false;
		for ( Plasmid p : targetRecipient.getPlasmidsHosted() )
			if ( (! plasmidCompatibilityMarkers.contains(p.getName())) || name.equals(p.getName()))
				return false;
		return true;    
	}
	/**
	 * \brief When plasmidID is given, returns String plasmidName
	 */
	public static String getPlasmidName (int plasmidID){
		return plasmidIDAndNameList.get(plasmidID);
	}
}
