package simulator.geometry.shape;

import java.io.Serializable;

import simulator.geometry.ContinuousVector;
import simulator.geometry.DiscreteVector;
import simulator.geometry.Domain;
import simulator.geometry.pointProcess.Edge;
import simulator.geometry.pointProcess.Site;
import utils.ExtraMath;
import utils.XMLParser;

public class Hemispherical implements IsShape, CanPointProcess, Serializable
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * A point on the cylinder axis
	 */
	private DiscreteVector _dPointCenterBase;
	
	/**
	 * 
	 */
	private DiscreteVector _dVectorToApex;
	
	/**
	 * A vector, orthogonal to the vector along the axis, used as a reference
	 * for the angle theta.
	 */
	private DiscreteVector _dVectorRadiusV;
	
	/**
	 * Another vector, orthogonal to the vector along the axis, used for
	 * converting back to Cartesian coordinates.
	 */
	private DiscreteVector _dVectorRadiusW;
	
	
	private ContinuousVector _cPointCenterBase;
	
	
	private ContinuousVector _cVectorToApex;
	
	/**
	 * A vector, orthogonal to the vector along the axis, used as a reference
	 * for the angle theta.
	 */
	private ContinuousVector _cVectorRadiusV;
	
	/**
	 * Another vector, orthogonal to the vector along the axis, used for
	 * converting back to Cartesian coordinates.
	 */
	private ContinuousVector _cVectorRadiusW;
	
	/**
	 * The radius of this cylinder. Equivalent to _cVectorRadius.norm().
	 */
	private Double _radius;
	
	/**
	 * Whether the inside of the hemisphere is the inside (true) or the outside
	 * (false) of the domain. 
	 */
	private Boolean _interiorMatchesDomain;
	
	/**
	 * 
	 * @param dPointCenterBase
	 * @param dVectorToApex
	 * @param interiorMatchesDomain
	 */
	public Hemispherical(DiscreteVector dPointCenterBase,
				DiscreteVector dVectorToApex, Boolean interiorMatchesDomain)
	{
		_dPointCenterBase = new DiscreteVector(dPointCenterBase);
		_dVectorToApex = new DiscreteVector(dVectorToApex);
		_interiorMatchesDomain = interiorMatchesDomain;
	}
	
	/**
	 * 
	 */
	public void readShape(XMLParser shapeRoot, Domain aDomain)
	{
		_dPointCenterBase = new DiscreteVector(shapeRoot.getParamParser("pointCenter"));
		_dVectorToApex = new DiscreteVector(shapeRoot.getParamParser("vectorAxis"));
		
		_dVectorRadiusV = new DiscreteVector();
		_dVectorRadiusW = new DiscreteVector();
		_dVectorToApex.orthoVector(_dVectorRadiusV, _dVectorRadiusW);
		
		Double res = aDomain.getResolution();
		_cPointCenterBase = new ContinuousVector(_dPointCenterBase, res);
		_cVectorToApex = new ContinuousVector(_dVectorToApex, res);
		_radius = _cVectorToApex.norm();
		_cVectorRadiusV = new ContinuousVector(_dVectorRadiusV, res);
		_cVectorRadiusV.normalizeVector(_radius);
		_cVectorRadiusW = new ContinuousVector(_dVectorRadiusW, res);
		_cVectorRadiusV.normalizeVector(_radius);
		
		_interiorMatchesDomain = shapeRoot.getParamBool("interiorMatchesDomain");
	}

	/**
	 * 
	 */
	public Boolean isOutside(ContinuousVector point)
	{
		ContinuousVector baseToPoint = baseToPoint(point);
		// Work out if the point is inside the hemisphere.
		Boolean isInsideHS = ( _cVectorToApex.cosAngle(baseToPoint) > 0 )
										&& ( baseToPoint.norm() <= _radius);
		// Use the exclusive or (Xor).
		return isInsideHS ^ _interiorMatchesDomain;
	}
	
	/**
	 * 
	 */
	public void orthoProj(ContinuousVector ccIn, ContinuousVector ccOut)
	{
		
		
	}
	
	/**
	 * 
	 */
	public ContinuousVector getOrthoProj(ContinuousVector ccIn)
	{
		
		return null;
	}
	
	/**
	 * 
	 */
	public Double getDistance(ContinuousVector point)
	{
		ContinuousVector baseToPoint = baseToPoint(point);
		Double cosAngle = _cVectorToApex.cosAngle(baseToPoint);
		if ( cosAngle >= 0 )
			return baseToPoint.norm() - _radius;
		/* If the point is below  the hemisphere, use the cosine rule to find
		 * the distance from the circle at the base. Note that, since we want
		 * cos(cosAngle - pi/2) we need to convert it first. 
		 */
		Double angle = Math.abs(Math.acos(cosAngle));
		angle -= Math.PI/2;
		cosAngle = Math.cos(angle);
		Double distSq = baseToPoint.prodScalar(baseToPoint);
		Double out = distSq + ExtraMath.sq(_radius);
		out  += 2 * Math.sqrt(distSq) * cosAngle;
		return Math.sqrt(out);
	}
	
	/**
	 * 
	 * TODO use convertToPolar(point, _cVectorRadiusV, _cVectorToApex)
	 * 
	 * @param point 
	 * @return 3 Doubles: (0) angle on the circular base, (1) angle with the
	 * vector to the apex, (2) distance from the centre.
	 */
	private Double[] convertToPolar(ContinuousVector point)
	{
		Double[] out = new Double[3];
		ContinuousVector baseToPoint = baseToPoint(point);
		
		out[1] = _cVectorToApex.cosAngle(baseToPoint);
		out[2] = baseToPoint.norm();
		
		// Should be equivalent to Planar.orthoProj() for the base of the HS
		ContinuousVector pointOnPlane = new ContinuousVector(_cVectorToApex);
		pointOnPlane.times(-out[1]*out[2]);
		pointOnPlane.add(baseToPoint);
		
		out[0] = _cVectorRadiusV.angle(pointOnPlane);
		out[1] = Math.acos(out[1]);
		return out;
	}
	
	private Double[] convertToPolar(ContinuousVector point,
					ContinuousVector radialVector, ContinuousVector apexVector)
	{
		Double[] out = new Double[3];
		ContinuousVector baseToPoint = baseToPoint(point);
		
		out[1] = apexVector.cosAngle(baseToPoint);
		out[2] = baseToPoint.norm();
		
		// Should be equivalent to Planar.orthoProj() for the base of the HS
		ContinuousVector pointOnPlane = new ContinuousVector(apexVector);
		pointOnPlane.times(-out[1]*out[2]);
		pointOnPlane.add(baseToPoint);
		
		out[0] = radialVector.angle(pointOnPlane);
		out[1] = Math.acos(out[1]);
		return out;
	}
	
	private void convertToCartesian(Double[] coords, ContinuousVector out,
							ContinuousVector center, ContinuousVector apexU,
							ContinuousVector radialV, ContinuousVector radialW)
	{
		out.set(center);
		
		ContinuousVector temp = new ContinuousVector(radialV);
		temp.times(coords[2] * Math.cos(coords[0]) * Math.sin(coords[1]));
		out.add(temp);
		
		temp.set(radialW);
		temp.times(coords[2] * Math.sin(coords[0]) * Math.sin(coords[1]));
		out.add(temp);
		
		temp.set(apexU);
		temp.times(coords[2] * Math.cos(coords[1]));
		out.add(temp);
	}
	
	private ContinuousVector convertToCartesian(Double[] coords)
	{
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(coords, out, _cPointCenterBase, _cVectorToApex,
											_cVectorRadiusV, _cVectorRadiusW);
		return out;
	}
	
	/**
	 * 
	 * Assumes both points on the surface!
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	public Double distance(ContinuousVector point1, ContinuousVector point2)
	{
		return _radius * baseToPoint(point1).angle(baseToPoint(point2));
	}
	
	/**
	 * 
	 * @param point Any point in Cartesian space.
	 * @return	The vector from _cPointCenterBase to this point.
	 */
	private ContinuousVector baseToPoint(ContinuousVector point)
	{
		ContinuousVector out = new ContinuousVector(point);
		out.subtract(_cPointCenterBase);
		return out;
	}
	
	private ContinuousVector midpoint(Site site1, Site site2)
	{
		ContinuousVector baseToS1 = baseToPoint(site1);
		ContinuousVector baseToS2 = baseToPoint(site2);
		ContinuousVector orthog = baseToS1.crossProduct(baseToS2);
		Double[] s2 = convertToPolar(site2, baseToS1, orthog);
		// TODO change multiplier for weighting
		s2[0] *= 0.5;
		
		ContinuousVector newW = baseToS1.crossProduct(orthog);
		newW.turnAround();
		newW.normalizeVector(_radius);
		
		ContinuousVector out = new ContinuousVector();
		convertToCartesian(s2, out, _cPointCenterBase, orthog, baseToS1, newW);
		return out;
	}
	
	/**
	 * 
	 */
	@Override
	public Edge bisect(Site site1, Site site2)
	{
		// If the sites are co-localised, there's no unique bisector.
		if ( site1.equals(site2) )
			return null;
		Edge out = new Edge();
		// Set the regions on either side of the Edge.
		out.region[0] = site1;
		out.region[1] = site2;
		// Already implicitly the case, but stated explicitly for clarity.
		out.endPoint[0] = null;
		out.endPoint[1] = null;
		// Find the angle each of these makes on the circular base
		Double[] s1 = convertToPolar(site1);
		Double[] s2 = convertToPolar(site2);
		//out.coefficient[0] = ((s2[0] + s1[0])/2) % (2 * Math.PI);
		
		return out;
	}
}
