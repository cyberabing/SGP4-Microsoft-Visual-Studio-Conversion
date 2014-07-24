#include "SpaceObject.h"

SpaceObject::SpaceObject(void){
	/* Default constructor, does nothing special. */
};

SpaceObject::SpaceObject(std::vector<std::string> threeLineElement, int COV, double bstarMultiplier, int wgsEllipsoid, double objectRadius, char opsMode){
	/* Create an object that represents something orbiting the Earth - a convenience interface for SGP4 propagator by D. Vallado.
	@param threeLineElement - a vector of strings that cointains the object name (line 0) and the first and second lines of the TLE file for a given object. Can have trailing end of line characters and must have the leading 1 and 2.
	@param COV - type of covariance to be used for the objects. If 1 an uncertainty sphere 1:1:1 will be used and only maximum collision probability will be computed.
		If 2 a number of three-line elements will be read for every object and covariance will be estimated for those using the method of V.P. Osweiler. It will be assumed that
		3LEs are in the files called SSC.txt where SSC if the catalogue number of every object and that the 3LEs are sorted oldest to last (the most recent one is last).
	@param bstarMultiplier - a factor by which the B* coefficient of a satellite will be multiplied by, used to simulate atmospheric density changes.
	@param wgsEllispoid - specifies which WGS Earth ellipsoid to use, possible options are 72, 84 and 721 (old WGS72).
	@param objectRadius - radius of the object in metres, used to compute the probability of collisions.
	@param opsMode - operation mode of the SGP4, either i (improved) or a (afspc).
	*/
	char OPSMODE = opsMode; // Improved operation mode of SGP4, should result in a smoother behaviour.
	earthEllipsoid = wgsEllipsoid;
	BstarMultiplier = bstarMultiplier;

	/* Initialise SGP4-specific attributes. */
	revnum=0; elnum=0;
	
	// These are defined in SGP4EXT, will need to uncomment them when not using it.
	//const double deg2rad = PI/180.0;         // Conversion factor. 0.0174532925199433
	//const double xpdotp = 1440.0/(2.0*PI);  // Minutes in a dayper a full revolution = 229.1831180523293. Unit: min/rad
	sgp4Sat = elsetrec(); // SGP4 satellite struct that is used for propagation using the SGP4.
	
	/* Define the desired EGS ellipsoid. */
	if(wgsEllipsoid==84){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==72){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==721){
		gravityModel = wgs72;
	};
	getgravconst( gravityModel, tumin, mu, Re, xke, j2, j3, j4, j3oj2 ); // Get values for these variables.

	/* Convert the TLE into the expected format. */
	char TLE1[130];	char TLE2[130]; // TLE lines - char arrays of the appropriate length that are expected by sgp4init.

	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		strcpy_s(TLE1, threeLineElement.at(1).c_str());
		strcpy_s(TLE2, threeLineElement.at(2).c_str());
	#else
		std::strcpy(TLE1, threeLineElement.at(1).c_str());
		std::strcpy(TLE2, threeLineElement.at(2).c_str());
	#endif
	
	/* Set the implied decimal points since doing a formated read
     * fixes for bad input data values (missing, ...). */
    for (int j = 10; j <= 15; j++){
        if (TLE1[j] == ' ')
			TLE1[j] = '_';
	};
	if (TLE1[44] != ' ')
		TLE1[43] = TLE1[44];
	TLE1[44] = '.';
	if (TLE1[7] == ' ')
		TLE1[7] = 'U';
    if (TLE1[9] == ' ')
		TLE1[9] = '.';
	for (int j = 45; j <= 49; j++){
		if (TLE1[j] == ' ')
			TLE1[j] = '0';
	};
    if (TLE1[51] == ' ')
		TLE1[51] = '0';
    if (TLE1[53] != ' ')
		TLE1[52] = TLE1[53];
	TLE1[53] = '.';
    TLE2[25] = '.';
    for (int j = 26; j <= 32; j++){
        if (TLE2[j] == ' ')
			TLE2[j] = '0';
	};
    if (TLE1[62] == ' ')
        TLE1[62] = '0';
    if (TLE1[68] == ' ')
        TLE1[68] = '0';

	/* Parse the TLE. */
	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		sscanf_s(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification, sizeof(char),intldesg, 11*sizeof(char),&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf_s(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#else
		sscanf(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification,intldesg,&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#endif
	
	std::stringstream strstream; // Record the NORAD ID as a string as well.
	strstream << sgp4Sat.satnum; strstream >> NORAD_ID;
		
	/* Iniitalise the SGP4 propagator variables and the propagator itself for this object. */
	// Find the TLE epoch in Julian days.
    int year = 2000 + sgp4Sat.epochyr; // N.B. this won't work for historic TLEs from before 2000.

	int epMon, epDay, epHr, epMin; double epSec; // TLE epoch components.
	days2mdhms(year, sgp4Sat.epochdays, epMon, epDay, epHr, epMin, epSec);
	jday(year, epMon, epDay, epHr, epMin, epSec, sgp4Sat.jdsatepoch); 

	// Find no, ndot, nddot.
    sgp4Sat.no   = sgp4Sat.no / xpdotp; // When using SGP4EXT this will already be in correct units, anyway multiply by: rad/min
    sgp4Sat.nddot= sgp4Sat.nddot * pow(10.0, nexp);
    sgp4Sat.bstar= sgp4Sat.bstar * pow(10.0, ibexp) * bstarMultiplier; // Multiply by the factor that allows variations in solar activity to be synthesised.
    
	// Convert to sgp4 units.
    sgp4Sat.a    = pow( sgp4Sat.no*tumin , (-2.0/3.0) );
    sgp4Sat.ndot = sgp4Sat.ndot  / (xpdotp*1440.0);  //* ? * minperday
    sgp4Sat.nddot= sgp4Sat.nddot / (xpdotp*1440.0*1440);
    
	// Find standard orbital elements.
    sgp4Sat.inclo = sgp4Sat.inclo  * deg2rad;
    sgp4Sat.nodeo = sgp4Sat.nodeo  * deg2rad;
    sgp4Sat.argpo = sgp4Sat.argpo  * deg2rad;
    sgp4Sat.mo    = sgp4Sat.mo     * deg2rad;

	// Iniitlaise the SGP4 for this TLE.
	sgp4init( gravityModel, OPSMODE, sgp4Sat.satnum, sgp4Sat.jdsatepoch-2433281.5, sgp4Sat.bstar, sgp4Sat.ecco, sgp4Sat.argpo, sgp4Sat.inclo, sgp4Sat.mo, sgp4Sat.no, sgp4Sat.nodeo, sgp4Sat);

	/* Call the propagator to get the initial state. */
	double r0[3]; double v0[3]; // Create arrays that SGP4 expects
    sgp4(gravityModel, sgp4Sat,  0.0, r0,  v0); // Third argument is time since epoch in minutes - 0.0 gives the initial position and velocity.
	currentPos.resize( sizeof(r0)/sizeof(r0[0]) ); currentVelo.resize( sizeof(v0)/sizeof(v0[0]) ); // Resize the vectors to be of appropriate length.
	currentPos = std::vector<double>(r0, r0+3); currentVelo = std::vector<double>(v0, v0+3); // Record the values of the initial position and velocity in the vectors.
	currentEpochJDAY = sgp4Sat.jdsatepoch;
	TLEepochJDAY = sgp4Sat.jdsatepoch; // In Julian Days.

	/* Get other values of interest, mianly orbital elements. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(r0, v0, mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	SemiMajorAxis = SMA; Eccentricity = ECC; Inclination = INCL; MeanAnomaly = MA; LongAscendingNode = LAN; ArgumentOfPerigee = ARGP;
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.

	hardBodyRadius = objectRadius/1000.; // Default value in km.

	if(COV==2){
		ReadThreeLineElements(); // Read three line elements from a file - need them to evaluate covariance at any epoch.
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
	}else if(COV==1){
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));

		// Initialise the diagonal entries of the covariance matrices to simulate the desired aspect ratio.
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0; FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0;
	}
};

void SpaceObject::PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* veloPtr, double JDAY, bool updateCurrentState){
	/* Propagate the satellite to the specified Julain day and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param JDAY - Julian day at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	*/
	double r[3]; double v[3]; // Create arrays that SGP4 expects
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0;
    sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);
	
	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length here as well - they are Cartesian components.
			currentPos.at(i) = r[i];
			currentVelo.at(i) = v[i];
		};
	};
};

void SpaceObject::Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second, bool updateCurrentState){
	double JDAY = currentEpochJDAY; // Initialise with the last new location.
	jday(year, month, day, hour, minute, second, JDAY); // Convert desired epoch to Julian days.
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0; // And minutes since TLE epoch.
	
	double r[3]; double v[3]; // Propagate using arrays that SGP4 expects
	sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);

	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length - they are Cartesian components.
			currentPos[i] = r[i];
			currentVelo[i] = v[i];
		};
	};
};

void SpaceObject::CalculatePerigeeRadius(void){
	/* Update the currentPerigeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
};

void SpaceObject::CalculateApogeeRadius(void){
	/* Update the currentApogeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.
};

void SpaceObject::SetHardBodyRadius(double bodyRadius){
	/* Set the hardBodyRadius attribute to the specified value in m. */
	hardBodyRadius = bodyRadius/1000.; // Convert to km.
};

void SpaceObject::SetCovarianceMatrixRTC(std::vector< std::vector<double> > CovarianceMatrixRTC){
	/* Set the FullCovarianceMatrixRTC attribute. It contains a 6x6 position and velocity covariance matrix in km and km/sec and is defined in
	radial - in-track - cross-track reference frame at the epoch of this object's TLE. Also initialise a set of TLEs stored in CovarianceTLEs
	that will be used to propagate the covariance to any epoch in a Monte Carlo fashion.
	@param CovarianceMatrixRTC - 6x6 position and velocity covariance matrix in km and km/sec defined in the radial - in-track - cross-track reference frame at the same epoch as this SpaceObject, i.e. TLEepochJDAY.
	*/
	FullCovarianceMatrixRTC = CovarianceMatrixRTC; // Set the full 6x6 covariance matrix.

	// Also set the 3x3 position covariance that will be used to compute the collision probability.
	PositionCovarianceMatrixRTC.at(0).at(0) = CovarianceMatrixRTC.at(0).at(0); PositionCovarianceMatrixRTC.at(0).at(1) = CovarianceMatrixRTC.at(0).at(1); PositionCovarianceMatrixRTC.at(0).at(2) = CovarianceMatrixRTC.at(0).at(2);
	PositionCovarianceMatrixRTC.at(1).at(0) = CovarianceMatrixRTC.at(1).at(0); PositionCovarianceMatrixRTC.at(1).at(1) = CovarianceMatrixRTC.at(1).at(1); PositionCovarianceMatrixRTC.at(1).at(2) = CovarianceMatrixRTC.at(1).at(2);
	PositionCovarianceMatrixRTC.at(2).at(0) = CovarianceMatrixRTC.at(2).at(0); PositionCovarianceMatrixRTC.at(2).at(1) = CovarianceMatrixRTC.at(2).at(1); PositionCovarianceMatrixRTC.at(2).at(2) = CovarianceMatrixRTC.at(2).at(2);

	//TODO: should spawn representative TLEs here to propagate the covariance in an MC fashion - save time and only do it once. Or maybe even spawn the TLEs when needed, otherwise leave them.
};

void SpaceObject::ReadThreeLineElements( void ){
	/* Read three-line elements from a file called NORAD_ID.txt and save them to CovarianceTLEs. Use these to estimate the covariance using the method by
	V.P. Osweiler in method ComputeCovarianceOSW. Assume that the 3LEs are in the order oldest to newest i.e. that the most-recent one is last. */
	std::string tempString = NORAD_ID; // Assume that the 3LEs from which to estimate the covariance are saved in files called NORAD_ID.txt. Make a copy of NORAD_ID attribute not to change it.
	const char* fileName = tempString.append(".txt").c_str();
	std::ifstream TLEfileStream(fileName, std::ifstream::in); // File streams that are necessary.

	CovarianceTLEs = std::vector<elsetrec>(); // TLEs that are representative of the covariance matrix and are used to propagate it in a Monte Carlo fashion.
	
	std::string TLEline; // Currently read TLE line.
	std::vector<std::string> currentTLE; currentTLE.resize(3); // Assembled TLE (second and third lines) and the object name (first line).

	int counterTLEs = 0; // Current line of the TLE, 0, 1 or 2.
	/* Read all the three-line elements from the file, those will be used to estimate and propagate the covariance matrix. */
	while( std::getline(TLEfileStream, TLEline) ){
		std::istringstream iss(TLEline);
		currentTLE.at(counterTLEs) = TLEline; // Add a new line.
		counterTLEs+=1;
		if(counterTLEs == 3){ // Read a whole TLE.
			counterTLEs = 0; // Re-set the counter.

			/* Convert the TLE into the expected format. */
			char TLE1[130];	char TLE2[130]; // TLE lines - char arrays of the appropriate length that are expected by sgp4init.
			elsetrec temporarySgp4Sat; // elsetrec for the curretnly-read TLE.

			#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
				strcpy_s(TLE1, currentTLE.at(1).c_str());
				strcpy_s(TLE2, currentTLE.at(2).c_str());
			#else
				std::strcpy(TLE1, threeLineElement.at(1).c_str());
				std::strcpy(TLE2, threeLineElement.at(2).c_str());
			#endif
	
			/* Set the implied decimal points since doing a formated read
			 * fixes for bad input data values (missing, ...). */
			for (int j = 10; j <= 15; j++){
				if (TLE1[j] == ' ')
					TLE1[j] = '_';
			};
			if (TLE1[44] != ' ')
				TLE1[43] = TLE1[44];
			TLE1[44] = '.';
			if (TLE1[7] == ' ')
				TLE1[7] = 'U';
			if (TLE1[9] == ' ')
				TLE1[9] = '.';
			for (int j = 45; j <= 49; j++){
				if (TLE1[j] == ' ')
					TLE1[j] = '0';
			};
			if (TLE1[51] == ' ')
				TLE1[51] = '0';
			if (TLE1[53] != ' ')
				TLE1[52] = TLE1[53];
			TLE1[53] = '.';
			TLE2[25] = '.';
			for (int j = 26; j <= 32; j++){
				if (TLE2[j] == ' ')
					TLE2[j] = '0';
			};
			if (TLE1[62] == ' ')
				TLE1[62] = '0';
			if (TLE1[68] == ' ')
				TLE1[68] = '0';

			/* Parse the TLE. */
			#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
				sscanf_s(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&temporarySgp4Sat.satnum,&classification, sizeof(char),intldesg, 11*sizeof(char),&temporarySgp4Sat.epochyr,&temporarySgp4Sat.epochdays,&temporarySgp4Sat.ndot,&temporarySgp4Sat.nddot,&nexp,&temporarySgp4Sat.bstar,&ibexp,&numb,&elnum);
				sscanf_s(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&temporarySgp4Sat.satnum,&temporarySgp4Sat.inclo,&temporarySgp4Sat.nodeo,&temporarySgp4Sat.ecco,&temporarySgp4Sat.argpo,&temporarySgp4Sat.mo,&temporarySgp4Sat.no,&revnum);
			#else
				sscanf(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&temporarySgp4Sat.satnum,&classification,intldesg,&temporarySgp4Sat.epochyr,&temporarySgp4Sat.epochdays,&temporarySgp4Sat.ndot,&temporarySgp4Sat.nddot,&nexp,&temporarySgp4Sat.bstar,&ibexp,&numb,&elnum);
				sscanf(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&temporarySgp4Sat.satnum,&temporarySgp4Sat.inclo,&temporarySgp4Sat.nodeo,&temporarySgp4Sat.ecco,&temporarySgp4Sat.argpo,&temporarySgp4Sat.mo,&temporarySgp4Sat.no,&revnum);
			#endif

			std::stringstream strstream; // Record the NORAD ID as a string as well.
			strstream << temporarySgp4Sat.satnum; strstream >> NORAD_ID;
		
			/* Iniitalise the SGP4 propagator variables and the propagator itself for this object. */
			// Find the TLE epoch in Julian days.
			int TLEyear = 2000 + temporarySgp4Sat.epochyr; // N.B. this won't work for historic TLEs from before 2000.

			int epMon, epDay, epHr, epMin; double epSec; // TLE epoch components.
			days2mdhms(TLEyear, temporarySgp4Sat.epochdays, epMon, epDay, epHr, epMin, epSec);
			jday(TLEyear, epMon, epDay, epHr, epMin, epSec, temporarySgp4Sat.jdsatepoch); 

			// Find no, ndot, nddot.
			temporarySgp4Sat.no   = temporarySgp4Sat.no / xpdotp; //* rad/min
			temporarySgp4Sat.nddot= temporarySgp4Sat.nddot * pow(10.0, nexp);
			temporarySgp4Sat.bstar= temporarySgp4Sat.bstar * pow(10.0, ibexp) * BstarMultiplier; // Multiply by the factor that allows variations in solar activity to be synthesised.
			
			// Convert to sgp4 units.
			temporarySgp4Sat.a    = pow( temporarySgp4Sat.no*tumin , (-2.0/3.0) );
			temporarySgp4Sat.ndot = temporarySgp4Sat.ndot  / (xpdotp*1440.0);  //* ? * minperday
			temporarySgp4Sat.nddot= temporarySgp4Sat.nddot / (xpdotp*1440.0*1440);
			
			// Find standard orbital elements.
			temporarySgp4Sat.inclo = temporarySgp4Sat.inclo  * deg2rad;
			temporarySgp4Sat.nodeo = temporarySgp4Sat.nodeo  * deg2rad;
			temporarySgp4Sat.argpo = temporarySgp4Sat.argpo  * deg2rad;
			temporarySgp4Sat.mo    = temporarySgp4Sat.mo     * deg2rad;

			// Iniitlaise the SGP4 for this TLE.
			sgp4init( gravityModel, OPSMODE, temporarySgp4Sat.satnum, temporarySgp4Sat.jdsatepoch-2433281.5, temporarySgp4Sat.bstar, temporarySgp4Sat.ecco, temporarySgp4Sat.argpo, temporarySgp4Sat.inclo, temporarySgp4Sat.mo, temporarySgp4Sat.no, temporarySgp4Sat.nodeo, temporarySgp4Sat);
			
			// Find the state vector of this TLE at the epoch of this SpaceObject.
			double r[3]; double v[3]; // Create arrays that SGP4 expects.
			double minutesSinceEpoch = (temporarySgp4Sat.jdsatepoch-TLEepochJDAY)*1440.0; // Minutes since epoch of this TLE to propagate it to the epoch of this SpaceObject.
			sgp4(gravityModel, temporarySgp4Sat, minutesSinceEpoch, r,  v);

			// Find classical orbital elements of this TLE at the epoch of this SpaceObject.
			double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
			rv2coe(r, v, mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );

			//TODO: filter out the outliers here.
			CovarianceTLEs.push_back( temporarySgp4Sat ); // Save this TLE.
		} // If read a full three-line element.
	} // While still reading the 3LE file.
};

std::vector< std::vector<double> > SpaceObject::InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr){
    /* Compute a rotation matrix from an inertial frame of reference (TEME, ECI or similar) to 
    Radial - Transverse - In-track (also called Radial - Along-track - Cross-ctrack) frame.
    @param position, velcity - dimensional position and velocity in the inertial frame. Must have the same units, assumed km and km/sec.
    @return - 3x3 rotation matrix from inertial to RTC frames - DOES NOT INCLUDE TRANSLATION BETWEEN ORIGINS.
	*/
	/* Base vectors of the RTC frame expressed in the inertial frame. */
    std::vector<double> radialUnit_inertial = vectorMultiplyByScalar( positionPtr, 1.0/vectorMagnitude(positionPtr) ); // Unit radial vector expressed in inertial frame.
	
    std::vector<double> crossUnit_inertial = std::vector<double>(3, 0.0); // Cross-track unit vector expressed in the inertial reference frame.
	crossProduct( &crossUnit_inertial, positionPtr, velocityPtr);
	unitVector( &crossUnit_inertial ); // Make unit.
    
	std::vector<double> transverseUnit_inertial = std::vector<double>(3, 0.0); // Transverse unit vector expressed in the inertial reference frame.
	crossProduct( &transverseUnit_inertial, &radialUnit_inertial, &crossUnit_inertial);
	unitVector( &transverseUnit_inertial ); // Make unit.

	/* Base vectors of the inertial frame expressed in the inertial frame. */
    std::vector<double> xAxisUnitVectorInertial = std::vector<double>(3, 0.0); xAxisUnitVectorInertial.at(0)=1.0;
    std::vector<double> yAxisUnitVectorInertial = std::vector<double>(3, 0.0); yAxisUnitVectorInertial.at(1)=1.0;
    std::vector<double> zAxisUnitVectorInertial = std::vector<double>(3, 0.0); zAxisUnitVectorInertial.at(2)=1.0;

	/* Formulate the rotation matrix from the inertial to RTC frame. */
    std::vector< std::vector<double> > rotationMatrix = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) ); // Full 3x3 rotation matrix from inertial to RTC coordinate frames.
	rotationMatrix.at(0).at(0) = dotProduct( &radialUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(0).at(1) = dotProduct( &radialUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(0).at(2) = dotProduct( &radialUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(1).at(0) = dotProduct( &transverseUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(1).at(1) = dotProduct( &transverseUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(1).at(2) = dotProduct( &transverseUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(2).at(0) = dotProduct( &crossUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(2).at(1) = dotProduct( &crossUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(2).at(2) = dotProduct( &crossUnit_inertial, &zAxisUnitVectorInertial );

	return rotationMatrix;
};

void SpaceObject::ComputeCovarianceOSW( double epochJDAY ){
	/* Compute full 6x6 covariance matrix (position and velocity) of Two-Line Elements at the desired epoch and save it in FullCovarianceMatrixRTC.
	Also save the 3x3 position covariance in PositionCovarianceMatrixRTC to use it later in true collision probability computations.
	This is done using the approach described by V. Osweiler 2006 i.e. by propagating all the TLEs to the desired epoch,
	computing residuals w.r.t. to the state vector corresponding to this SpaceObject at that epoch and estimating the covariance based on those residuals.
	Use TLEs saved in CovarianceTLEs to estimate the covariance (generated by @see ReadThreeLineELements).
	@param epochJDAY - Julian Day at which to compute the covariance.
	@return FullCovarianceMatrixRTC - 6x6 position (km) and velocity (km/sec) covariance matrix with variances on the diagonal and covariance coefficients off-diagonal.
		Positions are in rows 1 to 3, velocities 4 to 6. It is given in the velocity - normal - co-normal (cross-track) reference frame computed for this SpaceObject..
	@return PositionCovarianceMatrixRTC - 3x3 position (km) covariance matrix with variances on the diagonal and covariance coefficients off-diagonal.
		Positions are in rows 1 to 3. It is given in the velocity - normal - co-normal (cross-track) reference frame computed for this SpaceObject.
	*/

	/* Most-recent state and derived quantities. */
	std::vector<double> mostRecentPos = std::vector<double>(3, 0.0);
	std::vector<double> mostRecentVelo = std::vector<double>(3, 0.0);
	PropagateJDAY(&mostRecentPos, &mostRecentVelo, epochJDAY, false); // Get the position of the primary TLE (this SpaceObject) at the epoch of interest - it defines the local reference frame etc.

	std::vector< std::vector<double> > inertialToRTCrotation = InertialToRTC(&mostRecentPos, &mostRecentVelo); // Rotation matrix from the inertial frame to the RTC frame established by the state of this SpaceObject at the epoch of interest.

	/* Process all the TLEs that define the covariance matrix. */
	std::vector< std::vector<double> >* residuals = new std::vector< std::vector<double> >();
	residuals->resize(6); // Each entry of residuals is a vector with cartesian position and velocity components.S

	std::vector<double> tempPosResidual = std::vector<double>(3, 0.0); // 3 Cartesian position and 3 velocity components.
	std::vector<double> tempVeloResidual = std::vector<double>(3, 0.0);

	double meanX = 0.0, meanY = 0.0, meanZ = 0.0, meanVx=0.0, meanVy=0.0, meanVz=0.0; // Mean values of the residuals.
	int noResiduals = 0; // Number of the residuals read.

	/* Compute the residuals. */
	double r[3]; double v[3]; // Create arrays that SGP4 expects to find states of all the covariance-defining TLEs at the epoch of interest.
	for(std::vector<int>::size_type i=0; i<CovarianceTLEs.size()-1; i++){
		double minutesSinceEpoch = (epochJDAY-TLEepochJDAY)*1440.0; // Minutes since the prime TLE's epoch.
		sgp4(gravityModel, CovarianceTLEs.at(i), minutesSinceEpoch, r,  v); 
		for(std::vector<int>::size_type j=0; j<tempPosResidual.size(); j++){
			tempPosResidual.at(j) = r[j]-mostRecentPos.at(j);
			tempVeloResidual.at(j) = v[j]-mostRecentVelo.at(j);
		}

		tempPosResidual = vectorMultiplyByMatrix(&tempPosResidual, &inertialToRTCrotation); // Rotate the new residual to VNC established by the prime TLE.
		tempVeloResidual = vectorMultiplyByMatrix(&tempVeloResidual, &inertialToRTCrotation); 

		meanX += tempPosResidual.at(0); // Record all the residuals found to compute the average at the end.
		meanY += tempPosResidual.at(1);
		meanZ += tempPosResidual.at(2);
		meanVx += tempVeloResidual.at(0);
		meanVy += tempVeloResidual.at(1);
		meanVz += tempVeloResidual.at(2);

		for(std::vector<int>::size_type j=0; j<tempPosResidual.size(); j++){
			residuals->at(j).push_back( tempPosResidual.at(j) ); // Record this residual.
			residuals->at(j+3).push_back( tempVeloResidual.at(j) );
		}
		noResiduals += 1;
	}

	meanX = meanX/noResiduals; meanY = meanY/noResiduals; meanZ = meanZ/noResiduals;
	meanVx = meanVx/noResiduals; meanVy = meanVy/noResiduals; meanVz = meanVz/noResiduals;

	std::vector< std::vector<double> > residualsTranspose = transposeMatrix( residuals ); // Transpose the matrix.

	std::vector< std::vector<double> > posVeloCovariance = matrixMultiplyByMatrix(residuals, &residualsTranspose);
	FullCovarianceMatrixRTC = matrixMultiplyByScalar(&posVeloCovariance, 1.0/noResiduals); // Divide through by the number of residuals to get the completed covariance matrix.

	// Extract only the position components of the covariance matrix.
	PositionCovarianceMatrixRTC.at(0).at(0) = FullCovarianceMatrixRTC.at(0).at(0); PositionCovarianceMatrixRTC.at(0).at(1) = FullCovarianceMatrixRTC.at(0).at(1); PositionCovarianceMatrixRTC.at(0).at(2) = FullCovarianceMatrixRTC.at(0).at(2);
	PositionCovarianceMatrixRTC.at(1).at(0) = FullCovarianceMatrixRTC.at(1).at(0); PositionCovarianceMatrixRTC.at(1).at(1) = FullCovarianceMatrixRTC.at(1).at(1); PositionCovarianceMatrixRTC.at(1).at(2) = FullCovarianceMatrixRTC.at(1).at(2);
	PositionCovarianceMatrixRTC.at(2).at(0) = FullCovarianceMatrixRTC.at(2).at(0); PositionCovarianceMatrixRTC.at(2).at(1) = FullCovarianceMatrixRTC.at(2).at(1); PositionCovarianceMatrixRTC.at(2).at(2) = FullCovarianceMatrixRTC.at(2).at(2);
};

SpaceObject::~SpaceObject(void){
 /* Deconstructor, does nothing special.*/
};
