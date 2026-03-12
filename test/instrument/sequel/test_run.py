#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2026 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import datetime
import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PACKAGE_DIR = os.path.dirname(TEST_DIR)
sys.path.append(PACKAGE_DIR)

from anacore.instrument.sequel.run import (
    etreeToDict, getMetadataFromConsensusReadset, getRunInfo,
    getRunNameFromSequelHifi, Run
)


class TestConsensusReadset(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_file = os.path.join(self.tmp_folder, unique_id + ".consensusreadset.xml")
        with open(self.tmp_file, "w") as writer:
            writer.write("""<?xml version="1.0" encoding="utf-8"?>
<pbds:ConsensusReadSet CreatedAt="2026-02-28T10:14:58.598Z" MetaType="PacBio.DataSet.ConsensusReadSet" Name="lab1_Run324_smrt1_HIV_ENV-Cell1 (CCS)" Tags="ccs" TimeStampedName="lab1_Run324_smrt1_HIV_ENV-Cell1 (CCS)-260228_101458598" UniqueId="8dc59429-5166-4379-8833-095d24426646" Version="3.0.1" xmlns:pbbase="http://pacificbiosciences.com/PacBioBaseDataModel.xsd" xmlns:pbdm="http://pacificbiosciences.com/PacBioDataModel.xsd" xmlns:pbmeta="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd" xmlns:pbpn="http://pacificbiosciences.com/PacBioPartNumbers.xsd" xmlns:pbrk="http://pacificbiosciences.com/PacBioReagentKit.xsd" xmlns:pbsample="http://pacificbiosciences.com/PacBioSampleInfo.xsd" xmlns="http://pacificbiosciences.com/PacBioDatasets.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://pacificbiosciences.com/PacBioDatasets.xsd" xmlns:pbds="http://pacificbiosciences.com/PacBioDatasets.xsd">
	<pbbase:ExternalResources>
		<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.ConsensusReadFile.ConsensusReadBamFile" ResourceId="m54249Ue_260227_104832.hifi_reads.bam" TimeStampedName="pacbio_consensusreadfile_consensusreadbamfile-260228_101458599" UniqueId="0e3c20ab-a50f-465c-88f5-313883c12f03" Version="3.0.1">
			<pbbase:FileIndices>
				<pbbase:FileIndex CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.Index.PacBioIndex" ResourceId="m54249Ue_260227_104832.hifi_reads.bam.pbi" TimeStampedName="pacbio_index_pacbioindex-260228_101458599" UniqueId="8ab3c225-d0e4-4102-aec5-cd64330387cb" Version="3.0.1" />
			</pbbase:FileIndices>
			<pbbase:ExternalResources>
				<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.FileTypes.json" ResourceId="m54249Ue_260227_104832.zmw_metrics.json.gz" TimeStampedName="pacbio_filetypes_json-260228_101458599" UniqueId="e9b93cf7-b6be-4924-b089-1b31883012f1" Version="3.0.1" />
				<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.FileTypes.JsonReport" ResourceId="m54249Ue_260227_104832.ccs_reports.json" TimeStampedName="pacbio_filetypes_jsonreport-260228_101458599" UniqueId="23524104-0bb8-442d-a6e3-785cf6bea6e0" Version="3.0.1" />
				<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.FileTypes.log" ResourceId="m54249Ue_260227_104832.ccs.log" TimeStampedName="pacbio_filetypes_log-260228_101458599" UniqueId="e79bd0a7-d1e1-4f23-b892-b355422f85cc" Version="3.0.1" />
				<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.FileTypes.txt" ResourceId="m54249Ue_260227_104832.ccs_reports.txt" TimeStampedName="pacbio_filetypes_txt-260228_101458599" UniqueId="9130b927-632f-4657-94bd-77943fdea9a3" Version="3.0.1" />
				<pbbase:ExternalResource CreatedAt="2026-02-28T10:14:58.599Z" MetaType="PacBio.SubreadFile.ChipStatsFile" ResourceId="m54249Ue_260227_104832.sts.xml" TimeStampedName="pacbio_subreadfile_chipstatsfile-260228_101458599" UniqueId="e3bf3803-b841-4328-b2eb-33c26ff7e620" Version="3.0.1" />
			</pbbase:ExternalResources>
		</pbbase:ExternalResource>
	</pbbase:ExternalResources>
	<pbds:DataSetMetadata>
		<pbds:TotalLength>9615356986</pbds:TotalLength>
		<pbds:NumRecords>4895589</pbds:NumRecords>
		<Collections xmlns="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd">
			<CollectionMetadata Context="m54249Ue_260227_104832" CreatedAt="2026-02-27T10:01:20.3Z" InstrumentId="54249Ue" InstrumentName="Sequel II Instrument" MetaType="CollectionMetadata" ModifiedAt="0001-01-01T00:00:00" Status="Ready" TimeStampedName="54249Ue-CollectionMetadata-2026-34-27T10:34:10.207Z" UniqueId="53fe3f6e-9f43-443f-a6a4-727eedc8bef8">
				<pbmeta:MultiJobId>2539</pbmeta:MultiJobId>
				<pbmeta:ConsensusReadSetRef UniqueId="8dc59429-5166-4379-8833-095d24426646" />
				<pbmeta:RunDetails>
					<pbmeta:TimeStampedName>r54249Ue_20260227_103410</pbmeta:TimeStampedName>
					<pbmeta:Name>Run324_Smrt1_270226</pbmeta:Name>
					<pbmeta:CreatedBy>lab1</pbmeta:CreatedBy>
					<pbmeta:WhenCreated>2026-02-27T10:01:20.3Z</pbmeta:WhenCreated>
					<pbmeta:StartedBy>unknown</pbmeta:StartedBy>
					<pbmeta:WhenStarted>0001-01-01T00:00:00</pbmeta:WhenStarted>
				</pbmeta:RunDetails>
				<pbmeta:WellSample CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="lab1_Run324_smrt1_HIV_ENV">
					<pbmeta:WellName>A01</pbmeta:WellName>
					<pbmeta:Application>ampliconLessThan3kb</pbmeta:Application>
					<pbmeta:Concentration>0</pbmeta:Concentration>
					<pbmeta:OnPlateLoadingConcentration>110</pbmeta:OnPlateLoadingConcentration>
					<pbmeta:InsertSize>2324</pbmeta:InsertSize>
					<pbmeta:IsoSeq>false</pbmeta:IsoSeq>
					<pbmeta:IsCCS>true</pbmeta:IsCCS>
					<pbmeta:AutoAnalysisEnabled>false</pbmeta:AutoAnalysisEnabled>
					<pbmeta:SampleReuseEnabled>false</pbmeta:SampleReuseEnabled>
					<pbmeta:StageHotstartEnabled>false</pbmeta:StageHotstartEnabled>
					<pbmeta:SizeSelectionEnabled>false</pbmeta:SizeSelectionEnabled>
					<pbmeta:UseCount>0</pbmeta:UseCount>
					<pbmeta:DNAControlComplex>Sequel® II DNA Internal Control Complex 3.1</pbmeta:DNAControlComplex>
					<pbmeta:SampleSetupId>00e6205b-ee7d-7b2c-4a8d-f66be39c941d</pbmeta:SampleSetupId>
					<pbsample:BioSamples xmlns="http://pacificbiosciences.com/PacBioSampleInfo.xsd">
						<pbsample:BioSample Name="lab1_Run324_smrt1_HIV_ENV" />
					</pbsample:BioSamples>
				</pbmeta:WellSample>
				<pbmeta:Automation Name="Diffusion">
					<pbbase:AutomationParameters xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="ReuseSample" SimpleValue="False" ValueDataType="String" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="ImmobilizationTime" SimpleValue="120" ValueDataType="Double" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="ExtendFirst" SimpleValue="True" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="ExtensionTime" SimpleValue="30" ValueDataType="Double" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="MovieLength" SimpleValue="900" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="PCDinPlate" SimpleValue="True" ValueDataType="Boolean" />
						<pbbase:AutomationParameter Name="UsePredictiveLoadingBuffer" SimpleValue="True" ValueDataType="Boolean" />
						<pbbase:AutomationParameter Name="CollectionNumber" SimpleValue="0" ValueDataType="Int32" />
						<pbbase:AutomationParameter Name="InsertSize" SimpleValue="2324" ValueDataType="Int32" />
						<pbbase:AutomationParameter Name="BindingKitPartNumber" SimpleValue="102-194-200" ValueDataType="String" />
						<pbbase:AutomationParameter Name="PPAEstimatedDuration" SimpleValue="345" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="PPAMinDuration" SimpleValue="240.105" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="PPAMaxDuration" SimpleValue="449.895" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="PPATimeout" SimpleValue="1830" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="CCSOverheadDurationMinutesPerUmyGBase" SimpleValue="0.2" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="TipSearchMaxDuration" SimpleValue="576" ValueDataType="Int32" />
						<pbbase:AutomationParameter Name="SNRCut" SimpleValue="2" ValueDataType="Double" />
						<pbbase:AutomationParameter Name="HQRFMethod" SimpleValue="M1" ValueDataType="String" />
					</pbbase:AutomationParameters>
				</pbmeta:Automation>
				<pbmeta:CollectionNumber>0</pbmeta:CollectionNumber>
				<pbmeta:CellIndex>0</pbmeta:CellIndex>
				<pbmeta:SetNumber>0</pbmeta:SetNumber>
				<pbmeta:CellPac Barcode="DA310359" Description="Individual 4 Pack containing 4 SMRT®Cells each containing 7.76 million ZMWs" ExpirationDate="2027-03-28" LotNumber="420828" MovieTimeGrade="LR" Name="Sequel® II SMRT® Cell 8M (4/tray)" PartNumber="101-389-001" Version="1.0">
					<pbmeta:ChipLayout xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">Spider_1p0_NTO</pbmeta:ChipLayout>
				</pbmeta:CellPac>
				<pbmeta:ControlKit Barcode="Lxxxxx102249500123199" ChipType="8mChip" Description="Sequel II DNA Internal Control Complex 3.1 contains a fixed template of 2 kb length, bound to Sequel DNA Polymerase 2.1 for use as an internal sequencing control on the Sequel System. Reagent quantities provide spike-in controls for a minimum of 24 samples." ExpirationDate="2099-12-31" LotNumber="Lxxxxx" Name="Sequel® II DNA Internal Control Complex 3.1" PartNumber="102-249-500" Tags="Control Kit, CCK" Version="1.0">
					<pbmeta:CustomSequence xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">&gt;left_adapter\nTAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA\n&gt;right_adapter\nTAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA\n&gt;custom_sequence\nTGTCTAGGTCATCTCAACGTAGCTTTGACATATAACTTATCTAAAGTAATCCCTGCACACCTGTATGCATTATGCTTGCTATACACGCGGACACAGGCATATCATTTATTTTTTGCCATGTCCATTAATTGTTCAATAATTTTACTCACGGTATTTAATTTGATGTTGTGTTATATAGAATTGGAATTAAACTTATAAGGATGCTTAGACGTTGCATTATAAAAGTTTATGTACTAAGTATTTAAGACATTGGCATATGATTATAGCTTGACATTATTAAAAATTAATTAATTAAATCTCACACAATACTTATTCAAGACATTTTTACTAAGATAACCAAAGGAATGCGAACAAAATAATACTTAAAATATAAGACTTAGAAGTAATATGATCCAATAGTTACATATAGTACACTAAGTTCCTAAATTATATAACTTTAAAACAAAGTTACGAAATTTGGAAATAATTTTATTTAATCATATTTTCATAATAATGAAATACTGTTTATTTCAGTGGCGAAAAGAGATAATACGATTTTATAGTGATAGAATATCCTTGAAATATCTAAAGATAAAATTAGAAACTTTCTCTTTTCGCTGTAAAGCTATATGACTTAAAAATAACTTATACGCAAAGTATATTGCAGTGGAAACCCAAGAGTATAGTAGCCATGTAATCTCGGGTTCGAAACTACACGCCGCGCACGTAGTCAGATGGTCTGAAACTTGTCTGGGGCTGTTTGTTGACGGATGGAGACTTCACTAAGTGGCGTCAGGCGATGCGCACACACGGGACTCAATCCCGTAGCATGTTATGTGTCGTTCGAAACTCGTGCGTTCGAGATTTACGCCACATTGCCGGCTGGTCCAAGGACGTTATCTACCAGATGATACGGTCCAATTCGTAAGTTTGACTCACATAGTCGCGAACCGCGGAGCTGGAGAACAATAATTACCGGATGATTAGTTGACCATACGCACTATCATGCTCCGTGACTCAGTTTCCGCCATGGAGTTCTCACAGCCCCGTGTGTACCATAACTGCAGTAAGTAAGGACCTTGTTCGGAGGCCGACTCGTATTTCATATGATCTTAGTCTCGCCACCTTATCGCACGAATTGGGGGTGTCTTTTAGCCGACTCCGGCACGATCCGCCGGGAAGTTACTCGACCAGTTGCGGGACGCCCTAGTATGTTCGTATTACGTTCGATGCGTAAGCACCCCAGAGATTTTTGGCGGACGTTTCGGTAAATCATAGTAGAACCGGAGCGGTAAAGCTATTGATAACACGCAGGGACGAGCCAGTCGTCTAAGCTCCTCAGGGGTACCGTTCGCCGGACTACAGCCTGTCCCCGGCGGCCGCAACTGGGCTGCGATCCAGCCCCCGCTCCAAAAGGATGACTCGACCTTGCGCCTCGCGTACTCTGCTCTCGAGCTGTCTCCGTGGGCAATGCCGGCTCACGCTGTGGGGAACCCTGGACGCCCGGGCCGAGCCGACGTGGCCCCGCCCAGGCCTTTTCGTCGATCGCAGCTATGTACCCTGTGCTGGCCAGCGCTACTGCGCCGGCCATTAGCGGTGCGCTCTCGACTCGGCCCCAACGTAGACGGCGTCGCTGGCCGGATTCAAAGAAGTGAGCTACTACCATCGCGTGACGCCCTGCGGGCCTGAGTAACCGTGCACGAAGGACACCCCGTTCGTGGCGGGGGTTGCCTCCGCGACGGTCGCCAACGTTGGGGGTCGGTGCATTCAGGCGACGAGGGACCGCTGGTTTCCGGAGAGCGGCCTGTGCTCACACAGGTGCGGTCCATGGGGCCTGTGGATCCGGTTCTCCCACGCGTAGCGCCGGCGTTAGCATGGACGCTAAATAAGTATACGCCGGCAAAGGGAGTGTAGGCCGGCCCGAGGGCAATCGCGGTTACCGGGGTGGGGGAGCTCCCCGCACCAGCCTTGATGTGGTGTGCGAGCG</pbmeta:CustomSequence>
				</pbmeta:ControlKit>
				<pbmeta:TemplatePrepKit Barcode="Lxxxxx102141700123199" ChipType="8mChip" Description="The SMRTbell® Prep Kit contains reagent supplies to perform SMRTbell library preparations of primer-annealed SMRTbell libraries for insert sizes ranging from 500 bp to over 20 kb." ExpirationDate="2099-12-31" LotNumber="Lxxxxx" MaxInsertSize="20000" MinInsertSize="500" Name="SMRTbell® Prep Kit 3.0" PartNumber="102-141-700" Tags="Template Prep Kit, TPK" Version="1.0">
					<pbmeta:LeftAdaptorSequence xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT</pbmeta:LeftAdaptorSequence>
					<pbmeta:LeftPrimerSequence xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">aacggaggaggagga</pbmeta:LeftPrimerSequence>
					<pbmeta:RightAdaptorSequence xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT</pbmeta:RightAdaptorSequence>
					<pbmeta:RightPrimerSequence xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd">aacggaggaggagga</pbmeta:RightPrimerSequence>
				</pbmeta:TemplatePrepKit>
				<pbmeta:BindingKit Barcode="Lxxxxx102194200123199" ChipType="8mChip" Description="The Sequel II Binding Kit 3.1 contains reagent supplies to bind prepared DNA template libraries to the Sequel II Polymerase 2.1 in preparation for HiFi sequencing. The result is a DNA polymerase template complex ready for use on the Sequel II and Sequel IIe Systems. Sequel II Polymerase 2.1 should be used with Sequel II Sequencing Kit 2.0. Reagent quantities support at least 24 binding reactions, and at least 3 SMRT Cells 8M per reaction (72 cells total), depending on use case, sample size and concentration. Note: Sequel II Binding Kit 3.1 is not recommended for inserts larger than 3 kb" ExpirationDate="2099-12-31" LotNumber="Lxxxxx" Name="Sequel® II Binding Kit 3.1" PartNumber="102-194-200" Tags="Binding Kit, BDK" Version="1.0" />
				<pbmeta:SequencingKitPlate Barcode="141058101826100090526" ChipType="8mChip" Description="The DNA Sequencing Kit contains a sequencing reagent plate with chemistry for single molecule real-time sequencing on the PacBio Sequel® Dev with Dynamic Loading. Reagent quantities support 4 sequencing reactions to be used in conjunction with SMRT® Cell 4Pac(s).  (4 Cells max/Each Row supplies reagents for 1 Sequel SMRT Cell)" ExpirationDate="2026-09-05" LotNumber="141058" MaxCollections="4" Name="Sequel® II Sequencing Plate 2.0 (4 rxn)" NumOseTubes="0" PartNumber="101-826-100" SupportsDynamicLoading="true" Tags="Sequencing Kit, SQK" Version="2.0">
					<pbmeta:ReagentTubes Barcode="040384100619600083129" ExpirationDate="2029-08-31" LotNumber="040384" Name="Sequel® SMRT®Cell Oil" PartNumber="100-619-600" xmlns="http://pacificbiosciences.com/PacBioReagentKit.xsd" />
				</pbmeta:SequencingKitPlate>
				<pbmeta:Primary>
					<pbmeta:AutomationName>SequelAlpha</pbmeta:AutomationName>
					<pbmeta:ConfigFileName>SqlPoC_SubCrf_2C2A-t2.xml</pbmeta:ConfigFileName>
					<pbmeta:SequencingCondition>DefaultPrimarySequencingCondition</pbmeta:SequencingCondition>
					<pbmeta:CCSOptions>
						<pbmeta:ExecutionMode>OnInstrument</pbmeta:ExecutionMode>
						<pbmeta:IncludeKinetics>false</pbmeta:IncludeKinetics>
					</pbmeta:CCSOptions>
					<pbmeta:OutputOptions>
						<pbmeta:ResultsFolder>r54249Ue_20260227_103410/1_A01/</pbmeta:ResultsFolder>
						<pbmeta:CollectionPathUri>/home/pacbio/data_storage/r54249Ue_20260227_103410/1_A01/</pbmeta:CollectionPathUri>
						<pbmeta:CopyFiles>
							<pbmeta:CollectionFileCopy>Fasta</pbmeta:CollectionFileCopy>
							<pbmeta:CollectionFileCopy>Bam</pbmeta:CollectionFileCopy>
						</pbmeta:CopyFiles>
						<pbmeta:Readout>Bases_Without_QVs</pbmeta:Readout>
						<pbmeta:MetricsVerbosity>Minimal</pbmeta:MetricsVerbosity>
						<pbmeta:TransferResource>
							<pbmeta:Id>srs-pbi-collections</pbmeta:Id>
							<pbmeta:TransferScheme>SRS</pbmeta:TransferScheme>
							<pbmeta:Name>Transfer_srs_CHU</pbmeta:Name>
							<pbmeta:Description>Transfer scheme to use when running on instrument 54249</pbmeta:Description>
							<pbmeta:DestPath>/home/pacbio/data_storage</pbmeta:DestPath>
						</pbmeta:TransferResource>
					</pbmeta:OutputOptions>
				</pbmeta:Primary>
				<pbmeta:Secondary>
					<pbmeta:AutomationName>DefaultSecondaryAutomationName</pbmeta:AutomationName>
					<pbbase:AutomationParameters>
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="SourceId" SimpleValue="smrtlink-ui" ValueDataType="String" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="CCS" SimpleValue="True" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="CCSMode" SimpleValue="OnInstrument" ValueDataType="String" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="AllReads" SimpleValue="False" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="CpG" SimpleValue="False" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="AAV" SimpleValue="False" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="Kinetics" SimpleValue="False" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="Heteroduplex" SimpleValue="True" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="Demultiplex" SimpleValue="False" ValueDataType="Boolean" />
						<pbbase:AutomationParameter CreatedAt="2026-02-27T10:01:20.3Z" ModifiedAt="0001-01-01T00:00:00" Name="DemultiplexMode" SimpleValue="None" ValueDataType="String" />
					</pbbase:AutomationParameters>
					<pbmeta:CellCountInJob>0</pbmeta:CellCountInJob>
					<pbmeta:BarcodesFasta />
				</pbmeta:Secondary>
				<pbmeta:UserDefinedFields>
					<pbmeta:DataEntities Name=" LIMS_IMPORT " SimpleValue="DefaultUserDefinedFieldLIMS" ValueDataType="String" xmlns="http://pacificbiosciences.com/PacBioBaseDataModel.xsd" />
				</pbmeta:UserDefinedFields>
				<pbmeta:ComponentVersions>
					<pbmeta:VersionInfo Name="ics" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="iui" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="chemistry" Version="11.0.1.159204" />
					<pbmeta:VersionInfo Name="pa" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="paws" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="ppa" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="realtime" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="transfer" Version="11.0.1.162970" />
					<pbmeta:VersionInfo Name="smrtlink" Version="11.0.0.146107" />
					<pbmeta:VersionInfo Name="smrtimisc" Version="11.0.0.146107" />
					<pbmeta:VersionInfo Name="smrtinub" Version="11.0.0.146107" />
					<pbmeta:VersionInfo Name="smrtlink-analysisservices-gui" Version="11.0.0.146107" />
					<pbmeta:VersionInfo Name="smrttools" Version="11.0.0.146107" />
				</pbmeta:ComponentVersions>
			</CollectionMetadata>
		</Collections>
	</pbds:DataSetMetadata>
</pbds:ConsensusReadSet>""")

    def tearDown(self):
        if os.path.exists(self.tmp_file):
            os.remove(self.tmp_file)

    def testGetRunInfo(self):
        expected = {
            "InstrumentID": "54249Ue",
            "RunID": "m54249Ue_260227_104832",
            "RunName": "Run324_Smrt1_270226",
            "CreatedBy": "lab1",
            # StartedBy unknown
            "StartDate": datetime.datetime(2026, 2, 27, 10, 1, 20, tzinfo=datetime.timezone.utc),
            "EndDate": datetime.datetime(2026, 2, 28, 10, 14, 58, tzinfo=datetime.timezone.utc),
            "WellName": "A01",
            "Application": "ampliconLessThan3kb",
            "IsoSeq": False,
            "IsCCS": True
        }
        observed = getRunInfo(self.tmp_folder)
        self.assertEqual(expected, observed)

    def testGetRunNameFromSequelHifi(self):
        expected = "Run324_Smrt1_270226"
        observed = getRunNameFromSequelHifi(self.tmp_folder)
        self.assertEqual(expected, observed)


    def testGetMetadataFromConsensusReadset(self):
        expected = {'key': 'ConsensusReadSet', 'attributes': {'CreatedAt': '2026-02-28T10:14:58.598Z', 'MetaType': 'PacBio.DataSet.ConsensusReadSet', 'Name': 'lab1_Run324_smrt1_HIV_ENV-Cell1 (CCS)', 'Tags': 'ccs', 'TimeStampedName': 'lab1_Run324_smrt1_HIV_ENV-Cell1 (CCS)-260228_101458598', 'UniqueId': '8dc59429-5166-4379-8833-095d24426646', 'Version': '3.0.1', '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation': 'http://pacificbiosciences.com/PacBioDatasets.xsd'}, 'ExternalResources': {'key': 'ExternalResources', 'ExternalResource': {'key': 'ExternalResource', 'attributes': {'CreatedAt': '2026-02-28T10:14:58.599Z', 'MetaType': 'PacBio.ConsensusReadFile.ConsensusReadBamFile', 'ResourceId': 'm54249Ue_260227_104832.hifi_reads.bam', 'TimeStampedName': 'pacbio_consensusreadfile_consensusreadbamfile-260228_101458599', 'UniqueId': '0e3c20ab-a50f-465c-88f5-313883c12f03', 'Version': '3.0.1'}, 'FileIndices': {'key': 'FileIndices', 'FileIndex': {'key': 'FileIndex', 'attributes': {'CreatedAt': '2026-02-28T10:14:58.599Z', 'MetaType': 'PacBio.Index.PacBioIndex', 'ResourceId': 'm54249Ue_260227_104832.hifi_reads.bam.pbi', 'TimeStampedName': 'pacbio_index_pacbioindex-260228_101458599', 'UniqueId': '8ab3c225-d0e4-4102-aec5-cd64330387cb', 'Version': '3.0.1'}, 'val': None}}, 'ExternalResources': {'key': 'ExternalResources', 'ExternalResource': {'key': 'ExternalResource', 'attributes': {'CreatedAt': '2026-02-28T10:14:58.599Z', 'MetaType': 'PacBio.SubreadFile.ChipStatsFile', 'ResourceId': 'm54249Ue_260227_104832.sts.xml', 'TimeStampedName': 'pacbio_subreadfile_chipstatsfile-260228_101458599', 'UniqueId': 'e3bf3803-b841-4328-b2eb-33c26ff7e620', 'Version': '3.0.1'}, 'val': None}}}}, 'DataSetMetadata': {'key': 'DataSetMetadata', 'TotalLength': {'key': 'TotalLength', 'val': '9615356986'}, 'NumRecords': {'key': 'NumRecords', 'val': '4895589'}, 'Collections': {'key': 'Collections', 'CollectionMetadata': {'key': 'CollectionMetadata', 'attributes': {'Context': 'm54249Ue_260227_104832', 'CreatedAt': '2026-02-27T10:01:20.3Z', 'InstrumentId': '54249Ue', 'InstrumentName': 'Sequel II Instrument', 'MetaType': 'CollectionMetadata', 'Status': 'Ready', 'TimeStampedName': '54249Ue-CollectionMetadata-2026-34-27T10:34:10.207Z', 'UniqueId': '53fe3f6e-9f43-443f-a6a4-727eedc8bef8'}, 'MultiJobId': {'key': 'MultiJobId', 'val': '2539'}, 'ConsensusReadSetRef': {'key': 'ConsensusReadSetRef', 'attributes': {'UniqueId': '8dc59429-5166-4379-8833-095d24426646'}, 'val': None}, 'RunDetails': {'key': 'RunDetails', 'TimeStampedName': {'key': 'TimeStampedName', 'val': 'r54249Ue_20260227_103410'}, 'Name': {'key': 'Name', 'val': 'Run324_Smrt1_270226'}, 'CreatedBy': {'key': 'CreatedBy', 'val': 'lab1'}, 'WhenCreated': {'key': 'WhenCreated', 'val': '2026-02-27T10:01:20.3Z'}, 'StartedBy': {'key': 'StartedBy', 'val': None}, 'WhenStarted': {'key': 'WhenStarted', 'val': None}}, 'WellSample': {'key': 'WellSample', 'attributes': {'CreatedAt': '2026-02-27T10:01:20.3Z', 'Name': 'lab1_Run324_smrt1_HIV_ENV'}, 'WellName': {'key': 'WellName', 'val': 'A01'}, 'Application': {'key': 'Application', 'val': 'ampliconLessThan3kb'}, 'Concentration': {'key': 'Concentration', 'val': '0'}, 'OnPlateLoadingConcentration': {'key': 'OnPlateLoadingConcentration', 'val': '110'}, 'InsertSize': {'key': 'InsertSize', 'val': '2324'}, 'IsoSeq': {'key': 'IsoSeq', 'val': 'false'}, 'IsCCS': {'key': 'IsCCS', 'val': 'true'}, 'AutoAnalysisEnabled': {'key': 'AutoAnalysisEnabled', 'val': 'false'}, 'SampleReuseEnabled': {'key': 'SampleReuseEnabled', 'val': 'false'}, 'StageHotstartEnabled': {'key': 'StageHotstartEnabled', 'val': 'false'}, 'SizeSelectionEnabled': {'key': 'SizeSelectionEnabled', 'val': 'false'}, 'UseCount': {'key': 'UseCount', 'val': '0'}, 'DNAControlComplex': {'key': 'DNAControlComplex', 'val': 'Sequel® II DNA Internal Control Complex 3.1'}, 'SampleSetupId': {'key': 'SampleSetupId', 'val': '00e6205b-ee7d-7b2c-4a8d-f66be39c941d'}, 'BioSamples': {'key': 'BioSamples', 'BioSample': {'key': 'BioSample', 'attributes': {'Name': 'lab1_Run324_smrt1_HIV_ENV'}, 'val': None}}}, 'Automation': {'key': 'Automation', 'attributes': {'Name': 'Diffusion'}, 'AutomationParameters': {'key': 'AutomationParameters', 'AutomationParameter': {'key': 'AutomationParameter', 'attributes': {'Name': 'HQRFMethod', 'SimpleValue': 'M1', 'ValueDataType': 'String'}, 'val': None}}}, 'CollectionNumber': {'key': 'CollectionNumber', 'val': '0'}, 'CellIndex': {'key': 'CellIndex', 'val': '0'}, 'SetNumber': {'key': 'SetNumber', 'val': '0'}, 'CellPac': {'key': 'CellPac', 'attributes': {'Barcode': 'DA310359', 'Description': 'Individual 4 Pack containing 4 SMRT®Cells each containing 7.76 million ZMWs', 'ExpirationDate': '2027-03-28', 'LotNumber': '420828', 'MovieTimeGrade': 'LR', 'Name': 'Sequel® II SMRT® Cell 8M (4/tray)', 'PartNumber': '101-389-001', 'Version': '1.0'}, 'ChipLayout': {'key': 'ChipLayout', 'val': 'Spider_1p0_NTO'}}, 'ControlKit': {'key': 'ControlKit', 'attributes': {'Barcode': 'Lxxxxx102249500123199', 'ChipType': '8mChip', 'Description': 'Sequel II DNA Internal Control Complex 3.1 contains a fixed template of 2 kb length, bound to Sequel DNA Polymerase 2.1 for use as an internal sequencing control on the Sequel System. Reagent quantities provide spike-in controls for a minimum of 24 samples.', 'ExpirationDate': '2099-12-31', 'LotNumber': 'Lxxxxx', 'Name': 'Sequel® II DNA Internal Control Complex 3.1', 'PartNumber': '102-249-500', 'Tags': 'Control Kit, CCK', 'Version': '1.0'}, 'CustomSequence': {'key': 'CustomSequence', 'val': '>left_adapter\nTAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA\n>right_adapter\nTAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA\n>custom_sequence\nTGTCTAGGTCATCTCAACGTAGCTTTGACATATAACTTATCTAAAGTAATCCCTGCACACCTGTATGCATTATGCTTGCTATACACGCGGACACAGGCATATCATTTATTTTTTGCCATGTCCATTAATTGTTCAATAATTTTACTCACGGTATTTAATTTGATGTTGTGTTATATAGAATTGGAATTAAACTTATAAGGATGCTTAGACGTTGCATTATAAAAGTTTATGTACTAAGTATTTAAGACATTGGCATATGATTATAGCTTGACATTATTAAAAATTAATTAATTAAATCTCACACAATACTTATTCAAGACATTTTTACTAAGATAACCAAAGGAATGCGAACAAAATAATACTTAAAATATAAGACTTAGAAGTAATATGATCCAATAGTTACATATAGTACACTAAGTTCCTAAATTATATAACTTTAAAACAAAGTTACGAAATTTGGAAATAATTTTATTTAATCATATTTTCATAATAATGAAATACTGTTTATTTCAGTGGCGAAAAGAGATAATACGATTTTATAGTGATAGAATATCCTTGAAATATCTAAAGATAAAATTAGAAACTTTCTCTTTTCGCTGTAAAGCTATATGACTTAAAAATAACTTATACGCAAAGTATATTGCAGTGGAAACCCAAGAGTATAGTAGCCATGTAATCTCGGGTTCGAAACTACACGCCGCGCACGTAGTCAGATGGTCTGAAACTTGTCTGGGGCTGTTTGTTGACGGATGGAGACTTCACTAAGTGGCGTCAGGCGATGCGCACACACGGGACTCAATCCCGTAGCATGTTATGTGTCGTTCGAAACTCGTGCGTTCGAGATTTACGCCACATTGCCGGCTGGTCCAAGGACGTTATCTACCAGATGATACGGTCCAATTCGTAAGTTTGACTCACATAGTCGCGAACCGCGGAGCTGGAGAACAATAATTACCGGATGATTAGTTGACCATACGCACTATCATGCTCCGTGACTCAGTTTCCGCCATGGAGTTCTCACAGCCCCGTGTGTACCATAACTGCAGTAAGTAAGGACCTTGTTCGGAGGCCGACTCGTATTTCATATGATCTTAGTCTCGCCACCTTATCGCACGAATTGGGGGTGTCTTTTAGCCGACTCCGGCACGATCCGCCGGGAAGTTACTCGACCAGTTGCGGGACGCCCTAGTATGTTCGTATTACGTTCGATGCGTAAGCACCCCAGAGATTTTTGGCGGACGTTTCGGTAAATCATAGTAGAACCGGAGCGGTAAAGCTATTGATAACACGCAGGGACGAGCCAGTCGTCTAAGCTCCTCAGGGGTACCGTTCGCCGGACTACAGCCTGTCCCCGGCGGCCGCAACTGGGCTGCGATCCAGCCCCCGCTCCAAAAGGATGACTCGACCTTGCGCCTCGCGTACTCTGCTCTCGAGCTGTCTCCGTGGGCAATGCCGGCTCACGCTGTGGGGAACCCTGGACGCCCGGGCCGAGCCGACGTGGCCCCGCCCAGGCCTTTTCGTCGATCGCAGCTATGTACCCTGTGCTGGCCAGCGCTACTGCGCCGGCCATTAGCGGTGCGCTCTCGACTCGGCCCCAACGTAGACGGCGTCGCTGGCCGGATTCAAAGAAGTGAGCTACTACCATCGCGTGACGCCCTGCGGGCCTGAGTAACCGTGCACGAAGGACACCCCGTTCGTGGCGGGGGTTGCCTCCGCGACGGTCGCCAACGTTGGGGGTCGGTGCATTCAGGCGACGAGGGACCGCTGGTTTCCGGAGAGCGGCCTGTGCTCACACAGGTGCGGTCCATGGGGCCTGTGGATCCGGTTCTCCCACGCGTAGCGCCGGCGTTAGCATGGACGCTAAATAAGTATACGCCGGCAAAGGGAGTGTAGGCCGGCCCGAGGGCAATCGCGGTTACCGGGGTGGGGGAGCTCCCCGCACCAGCCTTGATGTGGTGTGCGAGCG'}}, 'TemplatePrepKit': {'key': 'TemplatePrepKit', 'attributes': {'Barcode': 'Lxxxxx102141700123199', 'ChipType': '8mChip', 'Description': 'The SMRTbell® Prep Kit contains reagent supplies to perform SMRTbell library preparations of primer-annealed SMRTbell libraries for insert sizes ranging from 500 bp to over 20 kb.', 'ExpirationDate': '2099-12-31', 'LotNumber': 'Lxxxxx', 'MaxInsertSize': '20000', 'MinInsertSize': '500', 'Name': 'SMRTbell® Prep Kit 3.0', 'PartNumber': '102-141-700', 'Tags': 'Template Prep Kit, TPK', 'Version': '1.0'}, 'LeftAdaptorSequence': {'key': 'LeftAdaptorSequence', 'val': 'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT'}, 'LeftPrimerSequence': {'key': 'LeftPrimerSequence', 'val': 'aacggaggaggagga'}, 'RightAdaptorSequence': {'key': 'RightAdaptorSequence', 'val': 'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT'}, 'RightPrimerSequence': {'key': 'RightPrimerSequence', 'val': 'aacggaggaggagga'}}, 'BindingKit': {'key': 'BindingKit', 'attributes': {'Barcode': 'Lxxxxx102194200123199', 'ChipType': '8mChip', 'Description': 'The Sequel II Binding Kit 3.1 contains reagent supplies to bind prepared DNA template libraries to the Sequel II Polymerase 2.1 in preparation for HiFi sequencing. The result is a DNA polymerase template complex ready for use on the Sequel II and Sequel IIe Systems. Sequel II Polymerase 2.1 should be used with Sequel II Sequencing Kit 2.0. Reagent quantities support at least 24 binding reactions, and at least 3 SMRT Cells 8M per reaction (72 cells total), depending on use case, sample size and concentration. Note: Sequel II Binding Kit 3.1 is not recommended for inserts larger than 3 kb', 'ExpirationDate': '2099-12-31', 'LotNumber': 'Lxxxxx', 'Name': 'Sequel® II Binding Kit 3.1', 'PartNumber': '102-194-200', 'Tags': 'Binding Kit, BDK', 'Version': '1.0'}, 'val': None}, 'SequencingKitPlate': {'key': 'SequencingKitPlate', 'attributes': {'Barcode': '141058101826100090526', 'ChipType': '8mChip', 'Description': 'The DNA Sequencing Kit contains a sequencing reagent plate with chemistry for single molecule real-time sequencing on the PacBio Sequel® Dev with Dynamic Loading. Reagent quantities support 4 sequencing reactions to be used in conjunction with SMRT® Cell 4Pac(s).  (4 Cells max/Each Row supplies reagents for 1 Sequel SMRT Cell)', 'ExpirationDate': '2026-09-05', 'LotNumber': '141058', 'MaxCollections': '4', 'Name': 'Sequel® II Sequencing Plate 2.0 (4 rxn)', 'NumOseTubes': '0', 'PartNumber': '101-826-100', 'SupportsDynamicLoading': 'true', 'Tags': 'Sequencing Kit, SQK', 'Version': '2.0'}, 'ReagentTubes': {'key': 'ReagentTubes', 'attributes': {'Barcode': '040384100619600083129', 'ExpirationDate': '2029-08-31', 'LotNumber': '040384', 'Name': 'Sequel® SMRT®Cell Oil', 'PartNumber': '100-619-600'}, 'val': None}}, 'Primary': {'key': 'Primary', 'AutomationName': {'key': 'AutomationName', 'val': 'SequelAlpha'}, 'ConfigFileName': {'key': 'ConfigFileName', 'val': 'SqlPoC_SubCrf_2C2A-t2.xml'}, 'SequencingCondition': {'key': 'SequencingCondition', 'val': 'DefaultPrimarySequencingCondition'}, 'CCSOptions': {'key': 'CCSOptions', 'ExecutionMode': {'key': 'ExecutionMode', 'val': 'OnInstrument'}, 'IncludeKinetics': {'key': 'IncludeKinetics', 'val': 'false'}}, 'OutputOptions': {'key': 'OutputOptions', 'ResultsFolder': {'key': 'ResultsFolder', 'val': 'r54249Ue_20260227_103410/1_A01/'}, 'CollectionPathUri': {'key': 'CollectionPathUri', 'val': '/home/pacbio/data_storage/r54249Ue_20260227_103410/1_A01/'}, 'CopyFiles': {'key': 'CopyFiles', 'CollectionFileCopy': {'key': 'CollectionFileCopy', 'val': 'Bam'}}, 'Readout': {'key': 'Readout', 'val': 'Bases_Without_QVs'}, 'MetricsVerbosity': {'key': 'MetricsVerbosity', 'val': 'Minimal'}, 'TransferResource': {'key': 'TransferResource', 'Id': {'key': 'Id', 'val': 'srs-pbi-collections'}, 'TransferScheme': {'key': 'TransferScheme', 'val': 'SRS'}, 'Name': {'key': 'Name', 'val': 'Transfer_srs_CHU'}, 'Description': {'key': 'Description', 'val': 'Transfer scheme to use when running on instrument 54249'}, 'DestPath': {'key': 'DestPath', 'val': '/home/pacbio/data_storage'}}}}, 'Secondary': {'key': 'Secondary', 'AutomationName': {'key': 'AutomationName', 'val': 'DefaultSecondaryAutomationName'}, 'AutomationParameters': {'key': 'AutomationParameters', 'AutomationParameter': {'key': 'AutomationParameter', 'attributes': {'CreatedAt': '2026-02-27T10:01:20.3Z', 'Name': 'DemultiplexMode', 'SimpleValue': 'None', 'ValueDataType': 'String'}, 'val': None}}, 'CellCountInJob': {'key': 'CellCountInJob', 'val': '0'}, 'BarcodesFasta': {'key': 'BarcodesFasta', 'val': None}}, 'UserDefinedFields': {'key': 'UserDefinedFields', 'DataEntities': {'key': 'DataEntities', 'attributes': {'Name': ' LIMS_IMPORT ', 'SimpleValue': 'DefaultUserDefinedFieldLIMS', 'ValueDataType': 'String'}, 'val': None}}, 'ComponentVersions': {'key': 'ComponentVersions', 'VersionInfo': {'key': 'VersionInfo', 'attributes': {'Name': 'smrttools', 'Version': '11.0.0.146107'}, 'val': None}}}}}}
        observed = getMetadataFromConsensusReadset(self.tmp_file)
        self.assertEqual(expected, observed)


class TestRun(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in_folder = os.path.join(tmp_folder, unique_id + "_run")
        os.makedirs(self.tmp_in_folder)

    def tearDown(self):
        if os.path.exists(self.tmp_in_folder):
            for filename in os.listdir(self.tmp_in_folder):
                os.remove(os.path.join(self.tmp_in_folder, filename))
            os.rmdir(self.tmp_in_folder)

    def testIsCopied(self):
        run = Run(self.tmp_in_folder)
        with open(os.path.join(self.tmp_in_folder, "m54249Ue_260227_104832.zmw_metrics.json.gz"), "w") as writer:
            writer.write("")
        self.assertFalse(run.isCopied())
        with open(os.path.join(self.tmp_in_folder, "m54249Ue_260227_104832.transferdone"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isCopied())

    def testIsEnded(self):
        run = Run(self.tmp_in_folder)
        with open(os.path.join(self.tmp_in_folder, "m54249Ue_260227_104832.zmw_metrics.json.gz"), "w") as writer:
            writer.write("")
        self.assertFalse(run.isEnded())
        with open(os.path.join(self.tmp_in_folder, "m54249Ue_260227_104832.transferdone"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isEnded())

    def testIsSequenced(self):
        run = Run(self.tmp_in_folder)
        self.assertFalse(run.isSequenced())
        with open(os.path.join(self.tmp_in_folder, "m54249Ue_260227_104832.zmw_metrics.json.gz"), "w") as writer:
            writer.write("")
        self.assertTrue(run.isSequenced())


if __name__ == "__main__":
    unittest.main()
