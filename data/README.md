README

This data file is published by the Movebank Data Repository (www.datarepository.movebank.org). As of the time of publication, a version of the published animal tracking data set can be viewed on Movebank (www.movebank.org) in the study "Navigation in Scopoli's shearwaters (data from Pollonara et al. 2015)". Individual attributes in the data files are defined below and in the Movebank Attribute Dictionary, available at www.movebank.org/node/2381.

This data package includes the following data files:
Navigation in Scopoli's shearwaters (data from Pollonara et al. 2015)-argos.csv
Navigation in Scopoli's shearwaters (data from Pollonara et al. 2015)-gps.csv
Navigation in Scopoli's shearwaters (data from Pollonara et al. 2015)-reference-data.csv

These data are described in the following written publication:
Pollonara E, Luschi P, Guilford T, Wikelski M, Bonadonna F, Gagliardo A (2015) Olfaction and topography, but not magnetic cues, control navigation in a pelagic seabird: displacements with shearwaters in the Mediterranean Sea. Scientific Reports. doi:10.1038/srep16486

Data package citation:
Pollonara E, Luschi P, Guilford T, Wikelski M, Bonadonna F, Gagliardo A (2015) Data from: Olfaction and topography, but not magnetic cues, control navigation in a pelagic seabird: displacements with shearwaters in the Mediterranean Sea. Movebank Data Repository. doi:10.5441/001/1.tc76g560

-----------

Terms of Use
This data file is licensed by the Creative Commons Zero (CC0 1.0) license. The intent of this license is to facilitate the re-use of works. The Creative Commons Zero license is a "no rights reserved" license that allows copyright holders to opt out of copyright protections automatically extended by copyright and other laws, thus placing works in the public domain with as little legal restriction as possible. However, works published with this license must still be appropriately cited following professional and ethical standards for academic citation.

We highly recommend that you contact the data creator if possible if you will be re-using or re-analyzing data in this file. Researchers will likely be interested in learning about new uses of their data, might also have important insights about how to properly analyze and interpret their data, and/or might have additional data they would be willing to contribute to your project. Feel free to contact us at support@movebank.org if you need assistance contacting data owners.

See here for the full description of this license
http://creativecommons.org/publicdomain/zero/1.0

-----------

Data Attributes
These definitions come from the Movebank Attribute Dictionary, available at www.movebank.org/node/2381.

animal ID: An individual identifier for the animal, provided by the data owner. This identifier can be a ring number, a name, the same as the associated tag ID, etc. If the data owner does not provide an Animal ID, an internal Movebank animal identifier may sometimes be shown.
	example: 91876A, Gary
	same as: individual-local-identifier

animal life stage: The age class or life stage of the animal at the beginning of the deployment. Can be years or months of age or terms such as "adult", "subadult" and "juvenile". Units should be defined in the values (e.g. "2 years").
	example: juvenile, adult
	units: Any units should be defined in the remarks.

animal reproductive condition: The reproductive condition of the animal at the beginning of the deployment.
	example: non-reproductive, pregnant

Argos altitude: Altitude used for location calculation, Argos DIAG format file (definition from Argos User's Manual 2011)
	example: 27
	units: meters

Argos best level: Best signal strength, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: -117
	units: decibels (dB)

Argos calculated frequency: Calculated frequency, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 401.6732709
	units: hertz (Hz)
	same as: Argos calcul freq

Argos latitude 1: Solution 1. platform latitude in degrees and thousandths of degrees, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 19.493
	units: decimal degrees, WGS84 reference system
	same as: Argos lat1

Argos latitude 2: Solution 2. platform latitude in degrees and thousandths of degrees, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 14.773
	units: decimal degrees, WGS84 reference system
	same as: Argos lat2

Argos location class: The location class retrieved from Argos. Allowed values are 0, 1, 2, 3, A, B, and Z (definition from Argos User's Manual 2011).
	example: A
	units: none
	same as: Argos lc

Argos longitude 1: Solution 1. platform longitude in degrees and thousandths of degrees, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 99.712
	units: decimal degrees, WGS84 reference system
	same as: Argos lon1

Argos longitude 2: Solution 2. platform longitude in degrees and thousandths of degrees, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 120.286
	units: decimal degrees, WGS84 reference system
	same as: Argos lon2

Argos Nmessages: The number of messages received [to calculate location], Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 8
	units: none
	same as: Argos nb mes

Argos Nmessages 120: The number of messages received by the satellite at a signal strength greater than -120 decibels, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 2
	units: none
	same as: Argos nb mes 120

Argos NOPC: The number of plausibility checks successful (from 0-4), Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 3
	units: none

Argos pass duration: Time elapsed between the first and last message received by the satellite, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 118
	units: seconds

Argos quality indicator: This indicator gives information on the transmitter in terms of two digits, X and Y. X is the first digit and indicates residual error on the frequency calculation; Y is the second digit and indicates transmitter oscillator frequency drift between two satellite passes. Values provided in Argos DIAG format files (definition from Argos User's Manual 2011). Allowed values are
X=0: No calculation of residual frequency error (fewer than four messages received)
X=1,2,3: Unsatisfactory convergence of calculation
X=4: Residual frequency error > 1.5 Hz
X=5: 0.15 Hz < residual frequency error < 1.5 Hz
X=6: Residual frequency error < 0.15 Hz
Y=0: No check on transmit frequency drift, as the two results are more than 12 hours apart.
Y=1: Frequency discrepancy > 400 Hz Probably due to transmit frequency discrepancy, change of oscillator, etc.
Y=2: Previous location is less than 1/2 hour old. Frequency discrepancy > 30 Hz, i.e. F/F (over 10 min) >2.5 E-8
Y=3: Frequency drift > 4 Hz/minute, i.e. F/F (10 min) > 1.10-7
Y=4: Frequency drift < 4 Hz/minute, i.e. F/F (10 min) < 1.10-7
Y=5: Frequency drift < 2 Hz/minute, i.e. F/F (10 min) < 5.10-8
Y=6: Frequency drift < 1 Hz/minute, i.e. F/F (10 min) < 2.5 . 10-8
Y=7: Frequency drift < 0.4 Hz/minute, i.e. F/F (10 min) < 1.10-8
Y=8: Frequency drift < 0.2 Hz/minute, i.e. F/F (10 min) < 5.10-9
	example: 68
	units: none
	same as: Argos Iq

Argos sensor 1: The value of the first Argos sensor, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 229
	units: specific to the sensor

Argos sensor 2: The value of the second Argos sensor, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 42
	units: specific to the sensor

Argos sensor 3: The value of the third Argos sensor, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 3
	units: specific to the sensor

Argos sensor 4: The value of the fourth Argos sensor, Argos DIAG format file (definition from Argos User's Manual 2011).
	example: 63
	units: specific to the sensor

attachment type: The way a tag is attached to an animal. Values are chosen from a controlled list:
	collar: The tag is attached by a collar around the animal's neck.
	glue: The tag is attached to the animal using glue.
	harness: The tag is attached to the animal using a harness.
	implant: The tag is placed under the skin of the an animal.
	tape: The tag is attached to the animal using tape.
	other: user specified

deploy off timestamp: The timestamp when the tag deployment ended.
	example: 2009-10-01 12:00:00.000
	format: yyyy-MM-dd HH:mm:ss.sss
	units: UTC (Coordinated Universal Time) or GPS time, which is a few leap seconds different from UTC
	same as: deploy off date

deploy on latitude: The geographic latitude of the location where the animal was released (intended primarily for instances in which the animal release and tag retrieval locations have higher accuracy than those derived from sensor data).
	example: 27.3516
	units: decimal degrees, WGS84 reference system

deploy on longitude: The geographic longitude of the location where the animal was released (intended primarily for instances in which the animal release and tag retrieval locations have higher accuracy than those derived from sensor data).
	example: -97.3321
	units: decimal degrees, WGS84 reference system

deploy on timestamp: The timestamp when the tag deployment started.
	example: 2008-08-30 18:00:00.000
	format: yyyy-MM-dd HH:mm:ss.sss
	units: UTC (Coordinated Universal Time) or GPS time, which is a few leap seconds different from UTC
	same as: deploy on date

deployment end comments: a description of the end of a tag deployment, such as cause of mortality or notes on the removal and/or failure of tag.
	example: Data transmission stopped after 108 days. Cause unknown.

deployment end type: A categorical classification of the tag deployment end. Values are chosen from a controlled list:
	captured: The tag remained on the animal but the animal was captured or confined.
	dead: The deployment ended with the death of the animal that was carrying the tag.
	equipment failure: The tag stopped working.
	fall off: The attachment of the tag to the animal failed, and it fell of accidentally.
	other
	released: The tag remained on the animal but the animal was released from captivity or confinement.
	removal: The tag was purposefully removed from the animal.
	unknown: The deployment ended by an unknown cause.

deployment ID: A unique identifier for the deployment of a tag on animal, provided by the data owner. If the data owner does not provide a Deployment ID, an internal Movebank deployment identifier may sometimes be shown.
	example: Jane-Tag42

duty cycle: Remarks associated with the duty cycle of a tag during the deployment, describing the times it is on/off and the frequency at which it transmits or records data.
	example: it turns off during the night
	units: Any units should be defined in the remarks.

event ID: An identifier for the set of information associated with each record or event in a data set. A unique event ID is assigned to every time-location or other time-measurement record in Movebank.
	example: 6340565
	units: none

ground speed: The estimated ground speed between consecutive locations.
	example: 7.22
	units: meters per second

latitude (decimal degree): The geographic longitude of a location along an animal track as estimated by the processed sensor data. Positive values are east of the Greenwich Meridian, negative values are west of it.
	example: -121.1761111
	units: decimal degrees, WGS84 reference system
	same as: location lat

life stage: The age class or life stage of the animal at the beginning of the deployment. Can be years or months of age or terms such as "adult", "subadult" and "juvenile". Units should be defined in the values (e.g. "2 years").
	example: juvenile, adult
	units: Any units should be defined in the remarks.

longitude (decimal degree): The geographic longitude of a location along an animal track as estimated by the processed sensor data. Positive values are east of the Greenwich Meridian, negative values are west of it.
	example: -121.1761111
	units: decimal degrees, WGS84 reference system
	same as: location long

manipulation comments: Additional comments about the way in which the animal was manipulated during the deployment. Use manipulation type to define the general type of manipulation.
	example: Relocated from breeding colony on Smithers Island to release location at 70.02E, 21.21S

manipulation type: The way in which the animal was manipulated during the deployment. Additional details about the manipulation can be provided using manipulation comments. Values are chosen from a controlled list:
	confined: The animal's movement was restricted to within a defined area.
	none: The animal received no treatment other than the tag attachment.
	relocated: The animal was released from a site other than the one at which it was captured.
	manipulated other: The animal was manipulated in some other way, such as a physiological manipulation.

manually marked outlier: An event marked manually as an outlier, typically using the Event Editor in Movebank. Allowed values are TRUE or FALSE.

ring ID: A number or color scheme for a band or ring attached to the animal.
	example: 26225
	units: none
	same as: animal ring id

sensor type: The type of sensor with which data were collected. Values are chosen from a controlled list:
	acceleration: The sensor collects acceleration data.
	accessory measurements: The sensor collects accessory measurements, such as battery voltage.
	Argos Doppler Shift: The sensor is using Argos Doppler shift for determining position.
	bird ring: The animal is identified by a ring that has a unique ID.
	GPS: The sensor uses GPS to find location and stores these.
	natural mark: The animal is identified by a natural marking.
	radio transmitter: The sensor is a classical radio transmitter.
	solar geolocator: The sensor collects light levels, which are used to determine position (for processed locations).
	solar geolocator raw: The sensor collects light levels, which are used to determine position (for raw light-level measurements).

sex: The sex of the biological individual(s) represented in the Occurrence. Values are from a controlled list:
	m: male
	f: female

study: The name of the study in Movebank in which data are stored.

study site: The name of the deployment site, for example a field station or colony.
	example: Pickerel Island North

tag ID: A unique identifier for the tag, provided by the data owner. If the data owner does not provide a tag ID, an internal Movebank tag identifier may sometimes be shown.
	example: 2342, ptt_4532
	same as: tag local identifier

tag manufacturer name: The company or person that produced the tag.
	example: Holohil
	same as: manufacturer

tag mass: The mass of the tag.
	example: 24
	units: grams

tag model: The model of the tag.
	example: T61
	same as: model

tag readout method: The way the data are received from the tag. Values are chosen from a controlled list:
	satellite: Data are transferred via satellite.
	phone network: Data are transferred via a phone network, such as GSM or AMPS.
	other wireless: Data are transferred via another form of wireless data transfer, such as a VHF radio transmitter/receiver.
	tag retrieval: The tag must be physically retrieved in order to obtain the data.

taxon: The scientific name of the species on which the tag was deployed, as defined by the Integrated Taxonomic Information System (ITIS, www.itis.gov). If the species name can not be provided, this should be the lowest level taxonomic rank that can be determined and that is used in the ITIS taxonomy. Additional information can be provided using the term taxon detail.
	example: Buteo swainsoni
	same as: species, animal taxon, individual taxon canonical name

timestamp: The date and time a sensor measurement was taken.
	example: 2008-08-14 18:31:00.000
	format: yyyy-MM-dd HH:mm:ss.sss
	units: UTC (Coordinated Universal Time) or GPS time, which is a few leap seconds different from UTC

visible: Determines whether an event is visible on the Movebank Search map. Values are calculated automatically, with FALSE indicating that the event has been marked as an outlier by manually marked outlier or algorithm marked outlier. Allowed values are TRUE or FALSE.

-----------

More Information
For more information about this repository, see www.movebank.org/node/15294, the FAQ at www.movebank.org/node/2220, or contact us at support@movebank.org.