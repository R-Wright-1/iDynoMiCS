package simulator.agent.zoo;

import idyno.SimTimer;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.LocatedAgent;
import simulator.agent.SpecialisedAgent;
import simulator.agent.Species;
import simulator.geometry.ContinuousVector;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief Bacterium class that can host a number of Plasmids of differing
 * species. 
 * 
 * <p>Since this extends BactEPS, instances of this class may also produce and
 * excrete EPS.</p>
 * 
 * <p>This class is an amalgamation of the <b>EpiBac</b> class, written by
 * Brian Merkey, and the <b>MultiEpiBac</b> class, written by Sonia 
 * Martins.</p>
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidBac extends BactEPS
{
	/**
	 * Plasmids hosted by this bacterium.
	 */
	private LinkedList<Plasmid> _plasmidHosted = new LinkedList<Plasmid>();

	/**
	 * All agents of type Bacterium to be collected once per timestep for use in chemostats only
	 * For efficiency and for preventing newly born agents to conjugate
	 * Only one instance of the fields below shared by all PlasmidBac
	 */
	private static LinkedList<Bacterium> _allBact = new LinkedList<Bacterium>();
	private static int _now;
	private static int _allCellsButMe;
	
	public static int _numTry;
	
	public static int _numTrans;
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidBac()
	{
		super();
		_speciesParam = new PlasmidBacParam();
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException
	{
		PlasmidBac out = (PlasmidBac) super.clone();
		out._plasmidHosted = new LinkedList<Plasmid>();
		Plasmid newPlasmid;
		for (Plasmid aPlasmid : _plasmidHosted )
		{
			newPlasmid = (Plasmid) aPlasmid.clone();
			out._plasmidHosted.add(newPlasmid);
		}
		return out;
	}
	
	/**
	 * Called during species creation to build the progenitor.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSimulator,
													XMLParser aSpeciesRoot)
	{
		boolean isCreatedByDivision = false;
		/*
		 * Initialisation of the BactEPS, and its superclasses.
		 */
		super.initFromProtocolFile(aSimulator, aSpeciesRoot);
		//Qian 10.2016: add following part. Otherwise, there will be no plasmid since the beginning of simulation.
		/*
		 * Create hosted plasmids.
		 */
		for (String aSpeciesName : aSpeciesRoot.getChildrenNames("plasmid"))
			initPlasmid(aSpeciesName, isCreatedByDivision);
		/*
		 * Genealogy and size management.
		 */
		init();
		/*
		 * Finally, grab all the plasmid species names for reporting.
		 */
		this.collectPlasmidSpeciesNames(aSimulator);
	}


	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData, boolean createdByDivision)
	{
		boolean isCreatedByDivision = false;
		/*
		 * First, grab all the plasmid species names for reporting.
		 */
		this.collectPlasmidSpeciesNames(aSim);
		/*
		 * Find the position to start at by using length and number of values
		 * read.
		 */
		int nValsRead = 3*getSpeciesParam().potentialPlasmids.size();
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * Read in info from the result file IN THE SAME ORDER AS IT WAS
		 * OUTPUT. HGT parameters:
		 */
		double r, d;
		int nCopy, spCounter = 0;
		for ( String plasmidName : getPotentialPlasmidNames() )
		{
			LogFile.writeLog("has been stepped");
			nCopy = Integer.parseInt(singleAgentData[iDataStart+3*spCounter]);
			if ( nCopy <= 0 )
				continue;
			r = Double.parseDouble(singleAgentData[iDataStart+3*spCounter+1]);
			d = Double.parseDouble(singleAgentData[iDataStart+3*spCounter+2]);
			Plasmid aPlasmid = this.initPlasmid(plasmidName, isCreatedByDivision);
			aPlasmid.setDetails(nCopy, r, d);
		}
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for ( int i = 0; i < iDataStart; i++ )
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData, createdByDivision);
	}
	
	@Override
	public PlasmidBac sendNewAgent() throws CloneNotSupportedException
	{
		PlasmidBac baby = (PlasmidBac) this.clone();
		updateSize();
		return baby;
	}
	
	@Override
	public void makeKid(boolean isCreatedByDivision) throws CloneNotSupportedException
	{
		/*
		 * Create the new instance and update the lineage.
		 */
		PlasmidBac baby = sendNewAgent();
		//LogFile.writeLog("before this._genealogy: " + this._genealogy.length() + " baby._genealogy: " + baby._genealogy.length());
		recordGenealogy(this, baby);
		//LogFile.writeLog("after this._genealogy: " + this._genealogy.length() + " baby._genealogy: " + baby._genealogy.length());
		
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * Compute movement to apply to both cells.
		 */
		setDivisionDirection(getInteractDistance(baby) / 2);
		baby._movement.subtract(_divisionDirection);
		_movement.add(_divisionDirection);
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth(isCreatedByDivision);
		/*
		 * Both daughter cells have an identical list of plasmids hosted.
		 * Loss-at-division could happen to either but not both, so first
		 * determine which daughter it could happen to, then determine if it
		 * does actually happen.
		 */
		for ( int i = 0; i < this._plasmidHosted.size(); i++ )
			if ( ExtraMath.getUniRandDbl() < 0.5 )
				this._plasmidHosted.get(i).applySegregation();
			else
				baby._plasmidHosted.get(i).applySegregation();
	}
	
	/**
	 * \brief Create a new PlasmidBac agent (who a priori is registered in at
	 * least one container).
	 * 
	 * <p>This agent is located on the relevant grid, and may host plasmids if
	 * stated in the protocol file.</p>
	 * 
	 * @param position Where to put this new agent.
	 * @param root XMLParser used in determining if the cell should contain a
	 * plasmid(s).
	 * @see {@link #createNewAgent(ContinuousVector)}
	 */
	public void createNewAgent(ContinuousVector position, XMLParser root, boolean isCreatedByDivision)
	{
		try 
		{
			// Get a clone of the progenitor.
			PlasmidBac baby = (PlasmidBac) sendNewAgent();
			baby.setFamily();
			baby.updateMass();
			
			/* If no mass defined, use the division radius to find the mass */
			// Note this should have been done already in initFromProtocolFile
			if ( this._totalMass == 0.0 )
			{
				guessMass();
				LogFile.writeLog("Warning: PlasmidBac.createNewAgent calling guessMass()");
			}
			// randomise its mass
			baby.randomiseMass();
			//System.out.println("RADIUS AT THIS POINT: "+this._totalRadius);
			baby.updateSize();
			
			this._myDivRadius = getDivRadius();
			baby._myDivRadius = getDivRadius();
			baby._myDeathRadius = getDeathRadius();
			
			baby.setLocation(position);
			baby.registerBirth(isCreatedByDivision);
			
			/*
			 * This is the part specific to PlasmidBac!
			 */
			
			Plasmid plasmid;
			for ( String plName : root.getChildrenNames("plasmid") )
			{
				plasmid = baby.initPlasmid(plName, isCreatedByDivision);
				plasmid.setDetails(1, SimTimer.getCurrentTime(), -Double.MAX_VALUE);
			}
		} 
		catch (CloneNotSupportedException e) 
		{
			LogFile.writeError(e, "PlasmidBac.createNewAgent()");
		}
	}
	
	/*************************************************************************
	 * BASIC GETTERS & SETTERS
	 ************************************************************************/
	
	@Override
	public PlasmidBacParam getSpeciesParam()
	{
		return (PlasmidBacParam) _speciesParam;
	}
	
	/**
	 * \brief Provides a list of all Plasmids hosted by this PlasmidBac host.
	 * 
	 * @return LinkedList of Plasmid objects hosted by this PlasmidBac.
	 */
	public LinkedList<Plasmid> getPlasmidsHosted()
	{
		return _plasmidHosted;
	}
	
	/*************************************************************************
	 * 
	 ************************************************************************/
	
	@Override
	public void internalStep()
	{
		// make this list when the first agent is called, regardless of whether it has a plasmid or not
		// this ensures the list does not contain any newborns
		// but only if chemostat
		if ( Simulator.isChemostat ) makeListAllBacteria();

		/*
		 * Check if any plasmids have illegal copy numbers.
		 */
		checkMissingPlasmid();
		LogFile.writeLog("PlasmidBac.internalStep is called");
		/*
		 * BactEPS internalStep methods.
		 */
		grow();
		updateSize();
		manageEPS();
		if ( willDivide() )
			divide();
		if ( willDie() )
			die(true);
		/*
		 * Now try conjugating.
		 */
		conjugate();
	}
	
	/**
	 * \brief Remove any plasmids whose presence is no longer legal.
	 */
	protected void checkMissingPlasmid()
	{
		int nIter = this._plasmidHosted.size();
		Plasmid aPlasmid;
		for(int i = 0; i < nIter; i++)
		{
			aPlasmid = this._plasmidHosted.removeFirst();
			if ( aPlasmid.getCopyNumber() <= 0 )
				this.killPlasmid(aPlasmid);
			else
				this._plasmidHosted.addLast(aPlasmid);				
		}
		/*
		 * If any plasmids have been lost, refresh the reactions encoded by
		 * remaining plasmids: we should only remove those reactions  uniquely
		 * provided by the plasmid that was lost.
		 */
		if ( this._plasmidHosted.size() < nIter )
			this.refreshPlasmidReactions();
	}
	
	/**
	 * \brief Initialise a Plasmid to be hosted by this PlasmidBac.
	 * 
	 * <p>Note that the new Plasmid will have default copy number for the
	 * species it belongs to.</p>
	 * 
	 * @param plasmidName Species name of the new Plasmid.
	 */
	private Plasmid initPlasmid(String plasmidName, boolean isCreatedByDivision)
	{
		Plasmid aPlasmid = null;
		try
		{
			aPlasmid = (Plasmid) 
							_species.getSpecies(plasmidName).sendNewAgent();
			this.welcomePlasmid(aPlasmid);
			aPlasmid.registerBirth(isCreatedByDivision);
			
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "PlasmidBac.initPlasmid("+plasmidName+")");
		}
		return aPlasmid;
	}
	
	/**
	 * \brief Receive a Plasmid into this host.
	 * 
	 * <p><b>[Rob 31July2015]</b> Removed updates to conjugation time, etc:
	 * This is now handled by the donor plasmid.</p>
	 * 
	 * @param aPlasmid Plasmid to be hosted by this PlasmidBac.
	 */
	public void welcomePlasmid(Plasmid aPlasmid)
	{
		this._plasmidHosted.add(aPlasmid);
		this.addPlasmidReactions(aPlasmid);
	}
	
	/**
	 * \brief Tell a Plasmid that it has been lost, and so should die.
	 * 
	 * <p><b>[Rob 31July2015]</b> Changed the part about losing reactions so
	 * that we now refresh all plasmid-encoded reactions. This avoids the
	 * possibility of a reaction being lost when another hosted plasmid still
	 * encodes for it.</p>  
	 * 
	 * @param plasmid Plasmid that was hosted by this PlasmidBac, but has now
	 * been lost.
	 */
	private void killPlasmid(Plasmid aPlasmid)
	{
		aPlasmid.die();
		LogFile.writeLog("Plasmid "+aPlasmid.sendName()+
											" lost from "+this.sendName());
	}
	
	/**
	 * \brief Add all the reactions encoded by a given Plasmid.
	 * 
	 * @param aPlasmid Plasmid whose reactions should be conferred to this
	 * PlasmidBac host.
	 */
	private void addPlasmidReactions(Plasmid aPlasmid)
	{
		for ( int reacIndex : aPlasmid.getReactionsEncoded() )
			this.addActiveReaction(allReactions[reacIndex], true);
	}
	
	/**
	 * \brief Lose all the reactions encoded by a given Plasmid.
	 * 
	 * <p><b>[Rob 31July2015]</b> Beware of losing reactions that are also
	 * encoded by other plasmids: don't assume that each encoded reaction is
	 * unique to a particular plasmid species!</p>
	 * 
	 * @param aPlasmid Plasmid whose reactions were conferred to this
	 * PlasmidBac host.
	 */
	private void losePlasmidReactions(Plasmid aPlasmid)
	{
		for ( int reacIndex : aPlasmid.getReactionsEncoded() )
			this.removeReaction(allReactions[reacIndex]);
	}
	
	/**
	 * \brief Refresh all the reactions conferred to this PlasmidBac host by
	 * its hosted Plasmids. 
	 */
	private void refreshPlasmidReactions()
	{
		for ( Plasmid aPlasmid : this._plasmidHosted )
			this.losePlasmidReactions(aPlasmid);
		for ( Plasmid aPlasmid : this._plasmidHosted )
			this.addPlasmidReactions(aPlasmid);
	}
	
	/**
	 * \brief Ask all your Plasmids to conjugate, if they are ready.
	 */
	protected void conjugate()
	{
		/*
		 * No point continuing if there are no plasmids in this host.
		 */
		if ( this._plasmidHosted.isEmpty() )
			return;
		/*
		 * Build a neighbourhood including only non-self Bacteria. The methods
		 * for this differ between chemostat and biofilm simulations:
		 * -	In the biofilm, all non-self Bacteria within reach of the
		 * 		Plasmid's pilus should be included.
		 * Need to do this for each plasmid as parameters may differ
		 * 
		 * [Rob 31July2016] No need to shuffle this list: 
		 * LocatedAgent.pickNeighbour() uses a randomly generated value to pick
		 * a random Bacterium.
		 */
		
		for ( Plasmid aPlasmid : this._plasmidHosted )
			if ( aPlasmid.isReadyToConjugate() )
			{
				if ( Simulator.isChemostat ) {
					/* calculate probability of collision of this donor with other PlasmidBac agents,
					multiply this probability by getScaledTone,
					then work out how many others to attempt to conjugate with and put these into potentials
					 */
					LinkedList<Bacterium> partners = screenAllPartners(aPlasmid);
					//was buildNbh(aPlasmid); but need to separate chemostat from biofilm stuff
					// screenAllPartners() adds to the _testTally and tryToSendPlasmid() subtracts
					// iterate through all partners
					// instead of calling seachConjugation we only do one line from searchConjugation (tryToSendPlasmid)
					for (Bacterium aTarget: partners) {
						aPlasmid.tryToSendPlasmid(aTarget);
					}
				} else { // biofilm branch
					HashMap<Bacterium, Double> 
					potentialRecps = new HashMap<Bacterium, Double>();
					potentialRecps = buildNbh(aPlasmid, aPlasmid.getPilusRange());
					// searchConjugation calls Plasmid.updateTestTallyScaleScanRate() then calls tryToSendPlasmid() while it canScan()
					this.searchConjugation(aPlasmid, potentialRecps);
				}
			}
	}
	
	/**
	 * \brief Create a list of all Bacterium agents the first time PlasmidBac.screenAllPartners() is called
	 * This method is only called in chemostat simulations, for biofilms we create
	 * a list of neighbours for each PlasmidBac individually with buildNbh()
	 * The list will contain self, so self needs to be skipped in screenAllPartners()
	 */
	void makeListAllBacteria() {
//		LogFile.writeLog("_now " +_now + "  SimTimer.getCurrentIter() " + SimTimer.getCurrentIter() );

		if (_now != SimTimer.getCurrentIter()) {
			_now = SimTimer.getCurrentIter();
			_allBact.clear();
			_allCellsButMe = 0;
			_numTry = 0;
			_numTrans = 0;

			/* Collision rate calculations should not be able to consider whether the potentials for conjugation
			 * are of the right class (PlasmidBac) or not, so include all bacteria, but do not include EPS particles
			 * because EPS particles only represent EPS, they are not proper agents
			 */
			for ( SpecialisedAgent aSA : _agentGrid.agentList ) {
				if (aSA instanceof Bacterium)
				{
					_allCellsButMe++;
					_allBact.add((Bacterium) aSA);
				}
			}
			_allCellsButMe--;
			LogFile.writeLog("_agentGrid.agentList.size() " + _agentGrid.agentList.size() );
		}
		else
			return;	
	}

	/**
	 * This function is only called for chemostat simulations where collisions are random	
	 * @param aPlasmid
	 * @return thoseToScreen (other Bacterium agents that are potential conjugation partners)
	 */
	protected LinkedList<Bacterium> screenAllPartners(Plasmid aPlasmid)
	{

		// Work out how many cells a given plasmid donor collides with per time unit
		Double tmpCollCoeff = aPlasmid.getSpeciesParam().collisionCoeff;
		
		// dt is the duration of the time step
		double dt = SimTimer.getCurrentTimeStep();
		double chemostatVol = _species.domain.length_X * _species.domain.length_Y * _species.domain.length_Z;
		// not a 'proper' probability as could be larger than one, we deal with that below
		double probToScreen = tmpCollCoeff * _allCellsButMe * dt / chemostatVol;
		// scale by growth tone of this donor
		probToScreen *= this.getScaledTone();
//		LogFile.writeLog("tmpCollCoeff: " + tmpCollCoeff + "  allCellsButMe: " + _allCellsButMe);
//		LogFile.writeLog("dt: " + dt + "  chemostatVol: " + chemostatVol + "probToScreen: " + probToScreen);

		// _testTally accumulates remainders from past time steps, so we add probToScreen
		// in tryToSendPlasmid(), _testTally is decremented every attempt to send a plasmid, so don't do it here
		aPlasmid._testTally += probToScreen;
		int numScreen = (int) Math.floor(aPlasmid._testTally);
		
		//int numScreen = (int) Math.floor(probToScreen + aPlasmid._testTally); // _testTally accumulates remainders from past time steps
		//aPlasmid._testTally = probToScreen + aPlasmid._testTally - numScreen; // the new remainder is saved in _testTally
		
		/*int numScreen = (int) Math.floor(probToScreen);
		double remainder = probToScreen - numScreen;
		double randyDbl = ExtraMath.getUniRandDbl(0.0, 1.0);
		if (randyDbl < remainder) {
			numScreen++;
		}
		*/
		
//		LogFile.writeLog("probToScreen: " + probToScreen + "  numScreen: " + numScreen);
		
		// randomly pick nScreen bacteria
		// Note that we use a random number generator that has a very long period so it is impossible to get duplicate random integers
		LinkedList<Bacterium> out = new LinkedList<Bacterium>();
		Bacterium thoseToScreen;
		int randyInt;
		int myIndex = _allBact.indexOf(this);
		for (int i = 0; i < numScreen; i++) {
			randyInt = ExtraMath.getUniRandInt(_allCellsButMe);
			while (randyInt == myIndex) // exclude self from thoseToScreen
					randyInt = ExtraMath.getUniRandInt(_allCellsButMe);
			thoseToScreen = _allBact.get(randyInt);
			out.add(thoseToScreen);
		}

		return out;
	}
	
	/**
	 * \brief This is used only for biofilms, whereas in chemostats, screenAllPartners() is used
	 * Add all non-self Bacteria within reach of this to a HashMap of
	 * potential recipients.
	 * 
	 * <p>The double values in the HashMap correspond to the Bacterium's 
	 * probability of being selected at random. If the species parameter
	 * <i>scaleScanProb</i> is false (default) then these probabilities are all
	 * equal, but if it is true then they are scaled by the distance from the
	 * host (this cell).</p>
	 * 
	 * <p>Parameter <b>nbhRadius</b> is typically the pilus length.</p>
	 * 
	 * @param nbhRadius double length (in um) of the maximum cell 
	 * surface-surface distance for another Bacterium to be considered a
	 * neighbor.
	 */
	public HashMap<Bacterium, Double> buildNbh(Plasmid aPlasmid, double nbhRadius)
	{
		HashMap<Bacterium, Double> out = new HashMap<Bacterium, Double>();
		/*
		 * nbhRadius gives the distance OUTSIDE the donor agent that touches a
		 * recipient agent, and so we need to subtract the radii from
		 * getDistance() or add the radii to nbhRadius.
		 * 
		 * getDistance(aLocAgent) gets distance between cell centres.
		 */
		double donorRadius = this.getRadius(false);
		/*
		 * Find all potential neighbors in the neighboring grid elements by updating _myNeighbors
		 */
		this.getPotentialShovers(nbhRadius + donorRadius);
		/*
		 * Now remove agents that are too far away in Euclidean space rather than the grid
		 */
		double distance;
		double recipRadius;
		double probVar = 1.0;
		double cumulativeProb = 0.0;
		for ( LocatedAgent recip : _myNeighbors )
		{
			/*
			 * First check that the potential recipient is not the current
			 * host, and that it is a Bacterium (or subclass, e.g. PlasmidBac)
			 */
			if ( recip == this || ! (recip instanceof Bacterium) )
				continue;
			/*
			 * Now filter by the Euclidean distance between cell surfaces.
			 */
			recipRadius = recip.getRadius(false);
			distance = getDistance(recip) - donorRadius - recipRadius;
			if ( distance > nbhRadius )
				continue;
			/*
			 * Finally, add the cell, together with a probability variable.
			 * By default, all potential recipients are treated equally.
			 * If scaleScanProb is true, scale the probability by the distance
			 * from the donor (reasoning is similar to the intensity of 
			 * sunlight as a function of distance from the Sun's surface). 
			 */
			if ( getSpeciesParam().scaleScanProb )
				probVar = ExtraMath.sq( donorRadius / (donorRadius+distance));
			out.put((Bacterium) recip, probVar);
			cumulativeProb += probVar;
		}
		this._myNeighbors.clear();
		/*
		 * Now scale all probabilities so that they sum to one.
		 */
		scaleProbabilities(out, cumulativeProb);
		return out;
	}
	
	/**
	 * \brief Scale all probability variables in a given HashMap by 1/sum.
	 * 
	 * <p>Accepting <b>sum</b> as input, rather than calculating it directly,
	 * is a shortcut but depends on an accurate value of <b>sum</b>. 
	 * 
	 * @param hm HashMap<Bacterium, Double> linking each Bacterium with a
	 * probability variable.
	 * @param sum The sum of all these probability variables.
	 */
	private void scaleProbabilities(HashMap<Bacterium, Double> hm, double sum)
	{
		final double scaler = 1.0/sum;
		hm.replaceAll((b, p) -> {return p*scaler;});
	}
	
	/**
	 * \brief Randomly select a Bacterium from the list of potential
	 * recipients generated in buildNbh() (used in biofilms only)
	 * 
	 * <p>Assumes the sum total of probability variables to be one.</p>
	 */
	protected Bacterium pickPotentialRecipient(HashMap<Bacterium, Double> hm)
	{
		double rand = ExtraMath.getUniRandDbl();
		double counter = 0.0;
		Iterator<Bacterium> it = hm.keySet().iterator();
		Bacterium bac;
		while( it.hasNext() )
		{
			bac = it.next();
			counter += hm.get(bac);
			if ( counter > rand )
				return bac;
		}
		return null;
	}
	
	/**
	 * \brief Search for partners and try to send them a plasmid. Only used in biofilms.
	 * 
	 * @param aPlasmid A Plasmid, hosted by this PlasmidBac, that should try
	 * to conjugate with neighboring bacteria.
	 */
	public void searchConjugation(Plasmid aPlasmid,
										HashMap<Bacterium, Double> potentials)
	{
		/*
		 * If there is nobody to conjugate with, there is nothing more to do.
		 */
		if ( potentials.isEmpty() )
			return;
		/*
		 * First scales the plasmid's scan rate from its host's growth tone.
		 * The plasmid will calculate the number of neighbours it can
		 * look at per timestep.
		 * Then add this to the _testTally (there may be some overflow from the previous timestep).
		 * 
		 * This is not used in chemostat simulations, which require different treatment of collisions
		 */
		aPlasmid.updateTestTallyScaleScanRate(this.getScaledTone());
		/*
		 * Find a recipient(s) and try to send them a plasmid.
		 */
		while ( aPlasmid.canScan() )
		{
			LogFile.writeLog("try to search conjugation");
			aPlasmid.tryToSendPlasmid(pickPotentialRecipient(potentials));
		}
	}
	
	/**
	 * \brief Growth tone as a linear interpolation between the two cutoffs
	 * specified in the protocol file. If no cutoffs specified, growth tone will be 1.0
	 * 
	 * <p>See Merkey <i>et al</i> (2011) p.5 for more details:
	 * <a href=
	 * "http://onlinelibrary.wiley.com/doi/10.1111/j.1462-2920.2011.02535.x/full#ss18"
	 * >link</a>.</p>
	 * 
	 * @return double value in the range [0.0, 1.0]
	 */
	public double getScaledTone()
	{
		// Default is no growth rate dependence of plasmid transfer rate
		// If low and high tonus cutoffs are not set in the protocol file
		// they are set to -Double.MAX_VALUE by default
		// In this case this function should return 1
		Double lowTonus = getSpeciesParam().lowTonusCutoff; 
		Double highTonus = getSpeciesParam().highTonusCutoff;
		Double theTonus = growthTone();
		Double scaledTone = 1.0; // default

		/*
		 * Too low, so return zero.
		 */
		if ( theTonus < lowTonus )
			scaledTone = 0.0;
		/*
		 * Middle case, so do linear interpolation.
		 */
		else if ( theTonus < highTonus )
			scaledTone = (theTonus-lowTonus) / (highTonus-lowTonus);
		/*
		 * If neither of these is called we have a high tonus,
		 * so just return 1.0 (same effect as no growth dependence).
		 */
		return scaledTone;
	}
	
	/**
	 * \brief Net growth rate as a fraction of the maximum rate.
	 * 
	 * <p>See Merkey <i>et al</i> (2011) p.5 for more details:
	 * <a href=
	 * "http://onlinelibrary.wiley.com/doi/10.1111/j.1462-2920.2011.02535.x/full#ss18"
	 * >link</a>.</p>
	 * 
	 * @return Net growth rate divided by maximum growth rate.
	 */
	public Double growthTone()
	{
		return _netGrowthRate / getSpeciesParam().maxGrowthRate;
	}
	
	/*************************************************************************
	 * REPORTING
	 ************************************************************************/
	
	/**
	 * \brief Using the Simulator species list, collect the names of all
	 * Plasmid species that could be hosted by this PlasmidBac species.
	 * 
	 * <p>This list will be used when writing output and reading in from
	 * output agent_State file.</p>
	 * 
	 * @param aSim The Simulator this is running in.
	 */
	private void collectPlasmidSpeciesNames(Simulator aSim)
	{
		LogFile.writeLog("This is "+aSim.speciesDic.get(this.speciesIndex));
		for ( Species aSpecies : aSim.speciesList )
		{
			if ( ! ( aSpecies.getProgenitor() instanceof Plasmid ) )
				continue;
			LogFile.writeLog("Looking at "+aSpecies.speciesName);
			if ( ! ((Plasmid) aSpecies.getProgenitor()).isHostCompatible(this) )
			{
				LogFile.writeLog("\tHost not compatible");
				continue;
			}
			LogFile.writeLog("\tAdded!");
			getSpeciesParam().addPotentialPlasmidName(aSpecies.speciesName);
		}
	}
	
	/**
	 * \brief Provides a list of all plasmid species names that this
	 * PlasmidBac species could host.
	 * 
	 * @return ArrayList<String> of plasmid species names.
	 */
	private ArrayList<String> getPotentialPlasmidNames()
	{
//		LogFile.writeLog("a list of plasmids" +this.getSpeciesParam().potentialPlasmids);
		return this.getSpeciesParam().potentialPlasmids;
	}
	
	/**
	 * \brief Update the header for report output.
	 * 
	 * <p>For every plasmid species that could be hosted by this PlasmidBac
	 * species, the copy number, last reception time, and last donation time
	 * will always be reported, even if this is zero for some cells in this
	 * species.</p>
	 * 
	 * @see simulator.agent.LocatedAgent#sendHeader()
	 */
	@Override
	public StringBuffer sendHeader()
	{
		StringBuffer header = super.sendHeader();
		header.append(",plasmidID,tEntry,numHT,numVT,plasGenealogy");
		return header;
	}
	
	/**
	 * \brief Creates an output string of information generated on this
	 * particular agent.
	 * 
	 * Used in creation of results files.
	 * Writes the data matching the header file.
	 * 
	 * @return	String containing results associated with this agent.
	 */
	@Override
	public StringBuffer writeOutput()
	{
		//LogFile.writeLog("start of PlasmidBac.writeOutput");
		StringBuffer tempString = super.writeOutput();
		// String plasName; replaced plasName with plasID because parsing strings in the XML output when almost all output are numbers is a pain
		int plasID = -1;
		double tEntry = -1.0;
		int numHT = -1;
		int numVT = -1;
		String plasGenealogy = "-1";
		//LogFile.writeLog("plasGenealogy before for loop: " + plasGenealogy);
		for ( Plasmid aPlasmid : _plasmidHosted )
		{
			plasID = aPlasmid.getPlasmidID();
			tEntry = aPlasmid.getBirthday();
			numHT = aPlasmid.getNumHT();
			numVT=  aPlasmid.getGeneration();
			//LogFile.writeLog("All 4 in for loop: " +plasID + "," + tEntry + "," + numHT + "," + numVT);
			//plasGenealogy = aPlasmid.getGenealogy().toString();
			//LogFile.writeLog("plasGenealogy in for loop: " + plasGenealogy);
		}
		//LogFile.writeLog("plasGenealogy: " + plasGenealogy);
		tempString.append(","+plasID + "," + tEntry + "," + numHT + "," + numVT);// + "," + plasGenealogy);
		return tempString;
	}
	
	/*************************************************************************
	 * POV-RAY
	 ************************************************************************/
	
	@Override
	public String getName()
	{
		StringBuffer out = new StringBuffer( _species.speciesName );
		//TODO
		return out.toString();
	}
	
	@Override
	public Color getColor()
	{
		// TODO
		PlasmidBacParam param = getSpeciesParam();
		/*
		 * Recipients have no plasmid.
		 * Transconjugant received the plasmid after birth.
		 * Donor received the plasmid before/at birth.
		if ( _plasmidHosted.isEmpty() )
			return param.rColor;
		else if ( (_plasmidHosted_numHT == 0  )
			return param.tColor;
		else
		 */
			return param.dColor;
	}
	
	@Override
	public void writePOVColorDefinition(FileWriter fr) throws IOException 
	{
		//TODO
	}
}
