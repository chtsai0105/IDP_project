CREATE DATABASE /*!32312 IF NOT EXISTS*/ `IDP` /*!40100 DEFAULT CHARACTER SET latin1 */;

USE `IDP`;

--
-- Table structure for table `Aspergillus_fumigatus`
--

DROP TABLE IF EXISTS `organism`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `organism` (
  organism_id int NOT NULL auto_increment,
  species varchar(512) NOT NULL,
  filename varchar(1024) NOT NULL,
  phylum varchar(128) NULL,
  taxonomy_id int NULL,
  PRIMARY KEY (organism_id),
  KEY idx_species (species),
  KEY idx_phylum (phylum)
);
/*!40101 SET character_set_client = @saved_cs_client */;

DROP TABLE IF EXISTS `idp_sequence`;
CREATE TABLE `idp_sequence` (
  idp_id int NOT NULL auto_increment,
  sequence_id varchar(128) NOT NULL,
  organism_id int NOT NULL,
  idp_percentage DECIMAL(4,2),
  PRIMARY KEY (idp_id),
  UNIQUE KEY idx_species_seq (organism_id,sequence_id),
  KEY idx_percent (idp_percentage),
  FOREIGN KEY (organism_id) REFERENCES organism(organism_id)
);
/*!40101 SET character_set_client = @saved_cs_client */;
DROP TABLE IF EXISTS `idr`;
CREATE TABLE `idr` (
  idr_id int NOT NULL auto_increment,
  idp_id int NOT NULL,
  idr_start int NOT NULL,
  idr_end   int NOT NULL,
  PRIMARY KEY (idr_id),
  FOREIGN KEY (idp_id) REFERENCES idp_sequence(idp_id)
);
