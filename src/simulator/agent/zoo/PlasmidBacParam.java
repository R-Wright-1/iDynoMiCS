package simulator.agent.zoo;

import java.awt.Color;
import java.util.ArrayList;

import simulator.Simulator;
import utils.LogFile;
import utils.UnitConverter;
import utils.XMLParser;

/**
 * \brief TODO
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidBacParam extends BactEPSParam
{
	/**
	 * The maximum growth rate this PlasmidBac can achieve. Parameter for
	 * growth dependency of Plasmid scan rate.
	 */
	public double maxGrowthRate = 1.0;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public double lowTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Parameter for growth dependency: value specifies a transition point. By
	 * default, this takes a very large, negative number so that it has no
	 * effect.
	 */
	public double highTonusCutoff = - Double.MAX_VALUE;
	
	/**
	 * Whether or not to scale scan probabilities by distance from the host
	 * (only applies in biofilm simulations).
	 */
	public boolean scaleScanProb = false;
	
	public String hostID;
	/**
	 * Colours for POV-Ray output.
	 * 
	 * TODO define these for multiple hosted plasmids.
	 */
	public Color dColor, tColor, rColor;
	
	/**
	 * A list of Plasmid species that could be hosted by PlasmidBacs of this
	 * species. Useful in writing reports for output.
	 */
	public ArrayList<String> potentialPlasmids;
	
	private Simulator simulator;
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidBacParam()
	{
		
	}
	
	/**
	 * 
	 */
	@Override
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		simulator = aSim;
		double tempDbl;
		Boolean tempBool;
		String tempStr = "";
		StringBuffer tempUnit = new StringBuffer("");
		
		/*
		 * Get the growth dependency parameters.
		 */
		tempDbl = aSpeciesRoot.getParamDbl("maxGrowthRate", tempUnit);
		maxGrowthRate = tempDbl*UnitConverter.time(tempUnit.toString());
		/*
		tempDbl = getSpeciesParameterDouble("maxGrowthRate",
											aSpeciesRoot, speciesDefaults);
		maxGrowthRate = Double.isFinite(tempDbl) ? tempDbl : maxGrowthRate;
		*/
		
		tempDbl = getSpeciesParameterDouble("lowTonusCutoff",
											aSpeciesRoot, speciesDefaults);
		lowTonusCutoff = Double.isFinite(tempDbl) ? tempDbl : lowTonusCutoff;
		
		
		tempDbl = getSpeciesParameterDouble("highTonusCutoff", aSpeciesRoot,
															speciesDefaults);
		highTonusCutoff = Double.isFinite(tempDbl)? tempDbl : highTonusCutoff;
		
		tempBool = getSpeciesParameterBool("scaleScanProb", aSpeciesRoot,
															speciesDefaults);
		scaleScanProb = (tempBool == null) ? scaleScanProb : tempBool;
		if ( scaleScanProb && Simulator.isChemostat )
		{
			LogFile.writeLogAlways("Cannot scale scan probabilities by"+
				"distance in the chemostat! Setting scaleScanProb to false");
			scaleScanProb = false;
		}
		
		tempStr = getSpeciesParameterString ("hostID", aSpeciesRoot, speciesDefaults);
		tempStr = (tempStr==null) ? "TheOneAndOnly" : tempStr;
		hostID = tempStr;
						
		/*
		 * Get the colours for POV-Ray output.
		 */
		tempStr = getSpeciesParameterString("donorColor",
											aSpeciesRoot, speciesDefaults);
		tempStr = ( tempStr == null ) ? "white" : tempStr;
		dColor = utils.UnitConverter.getColor(tempStr);
		
		
		tempStr = getSpeciesParameterString("transconjugantColor",
											aSpeciesRoot, speciesDefaults);
		tempStr = ( tempStr == null ) ? "white" : tempStr;
		tColor = utils.UnitConverter.getColor(tempStr);
		
		
		tempStr = getSpeciesParameterString("recipientColor",
											aSpeciesRoot, speciesDefaults);
		tempStr = ( tempStr == null ) ? "white" : tempStr;
		rColor = utils.UnitConverter.getColor(tempStr);
		
		potentialPlasmids = new ArrayList<String>();
	}
	
	/**
	 * \brief Add the given name to the list of Plasmid species that could be
	 * hosted by PlasmidBacs of this species.
	 * 
	 * <p>Checks that the given name is not already on the list.</p>
	 * 
	 * @param name Plasmid species name.
	 */
	public void addPotentialPlasmidName(Plasmid aPlasmid)
	{
		if ( ! potentialPlasmids.contains(aPlasmid.getName()) ){
			potentialPlasmids.add(aPlasmid.getName());
			simulator.addPlasmidToGlobalList(aPlasmid);
		}
		/* for debugging
		for (int i =0;i<potentialPlasmids.size();i++)
		{
			LogFile.writeLog("potential plasmids: " + potentialPlasmids.get(i));
		}
		*/
	}
}
