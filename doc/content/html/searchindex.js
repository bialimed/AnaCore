Search.setIndex({docnames:["anacore","index"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.viewcode":1,sphinx:56},filenames:["anacore.rst","index.rst"],objects:{"anacore.STARLog":{STARLog:[0,1,1,""]},"anacore.abstractFile":{AbstractFile:[0,1,1,""],isEmpty:[0,3,1,""],isGzip:[0,3,1,""]},"anacore.abstractFile.AbstractFile":{close:[0,2,1,""],closed:[0,2,1,""],isRecordLine:[0,2,1,""],read:[0,2,1,""]},"anacore.annotVcf":{AnnotVCFIO:[0,1,1,""],VEPVCFIO:[0,1,1,""]},"anacore.annotVcf.AnnotVCFIO":{copyHeader:[0,2,1,""],recToVCFLine:[0,2,1,""],writeHeader:[0,2,1,""]},"anacore.bed":{BEDIO:[0,1,1,""],BEDRecord:[0,1,1,""],getAreas:[0,3,1,""],getAreasByChr:[0,3,1,""],getSortedAreasByChr:[0,3,1,""]},"anacore.bed.BEDIO":{BEDRecordToBEDLine:[0,2,1,""],getMaxNbCol:[0,2,1,""],isRecordLine:[0,2,1,""],isValid:[0,2,1,""],write:[0,2,1,""]},"anacore.bed.BEDRecord":{recFromRegion:[0,2,1,""]},"anacore.filters":{Filter:[0,1,1,""],FiltersCombiner:[0,1,1,""],filtersFromDict:[0,3,1,""]},"anacore.filters.Filter":{eval:[0,2,1,""],fromDict:[0,2,1,""],getRecordValue:[0,2,1,""],setAggregator:[0,2,1,""],setFct:[0,2,1,""],setGetter:[0,2,1,""],setOperator:[0,2,1,""],toDict:[0,2,1,""]},"anacore.filters.FiltersCombiner":{eval:[0,2,1,""]},"anacore.fusion":{ArribaIO:[0,1,1,""],BreakendVCFIO:[0,1,1,""],FusionCatcherIO:[0,1,1,""],FusionFileReader:[0,1,1,""],STARFusionIO:[0,1,1,""],getAltFromCoord:[0,3,1,""],getBNDInterval:[0,3,1,""],getCoordDictFromCoordStr:[0,3,1,""],getCoordStr:[0,3,1,""],getStrand:[0,3,1,""]},"anacore.fusion.ArribaIO":{isValid:[0,2,1,""],setVCFHeader:[0,2,1,""],titles:[0,4,1,""],write:[0,2,1,""]},"anacore.fusion.BreakendVCFIO":{close:[0,2,1,""],get:[0,2,1,""],isRecordLine:[0,2,1,""],isValid:[0,2,1,""],loadIndex:[0,2,1,""],write:[0,2,1,""]},"anacore.fusion.FusionCatcherIO":{isValid:[0,2,1,""],setVCFHeader:[0,2,1,""],titles:[0,4,1,""],write:[0,2,1,""]},"anacore.fusion.FusionFileReader":{factory:[0,2,1,""]},"anacore.fusion.STARFusionIO":{isValid:[0,2,1,""],setVCFHeader:[0,2,1,""],titles:[0,4,1,""],write:[0,2,1,""]},"anacore.genomicRegion":{CDS:[0,1,1,""],Exon:[0,1,1,""],Gene:[0,1,1,""],Intron:[0,1,1,""],Protein:[0,1,1,""],Transcript:[0,1,1,""]},"anacore.genomicRegion.Protein":{aaLength:[0,2,1,""],contains:[0,2,1,""],getCDSFromTranscript:[0,2,1,""],getNtPosFromRefPos:[0,2,1,""],getNtPosFromRegionPos:[0,2,1,""],getPosOnRef:[0,2,1,""],getPosOnRegion:[0,2,1,""],getSubFromRegionPos:[0,2,1,""],hasOverlap:[0,2,1,""],length:[0,2,1,""],setTranscript:[0,2,1,""]},"anacore.genomicRegion.Transcript":{addProtein:[0,2,1,""],getPosOnRef:[0,2,1,""],getPosOnRegion:[0,2,1,""],getSubFromRefPos:[0,2,1,""],getSubFromRegionPos:[0,2,1,""],length:[0,2,1,""],setProteins:[0,2,1,""]},"anacore.gff":{GFF3IO:[0,1,1,""],GFF3Record:[0,1,1,""]},"anacore.gff.GFF3IO":{decodedValue:[0,2,1,""],encodedValue:[0,2,1,""],isRecordLine:[0,2,1,""],write:[0,2,1,""]},"anacore.gff.GFF3Record":{addToAnnot:[0,2,1,""],toGff:[0,2,1,""]},"anacore.gtf":{GTFIO:[0,1,1,""],loadModel:[0,3,1,""]},"anacore.gtf.GTFIO":{isRecordLine:[0,2,1,""],recordToLine:[0,2,1,""],write:[0,2,1,""]},"anacore.hgvs":{Accession:[0,1,1,""],HGVS:[0,1,1,""],MutalyzerBatch:[0,1,1,""],RunMutalyzerDescription:[0,1,1,""],RunMutalyzerLegend:[0,1,1,""]},"anacore.hgvs.HGVS":{fromStr:[0,2,1,""],isPredicted:[0,2,1,""]},"anacore.hgvs.MutalyzerBatch":{getExecTime:[0,2,1,""],getParsedResponse:[0,2,1,""],getRequestURL:[0,2,1,""],getResponse:[0,2,1,""],submit:[0,2,1,""],syncRequest:[0,2,1,""],updateStatus:[0,2,1,""]},"anacore.hgvs.RunMutalyzerDescription":{getByAccession:[0,2,1,""]},"anacore.hgvs.RunMutalyzerLegend":{getIdByName:[0,2,1,""],getProtBytr:[0,2,1,""]},"anacore.illumina":{ADSSampleSheetIO:[0,1,1,""],CompletedJobInfo:[0,1,1,""],RTAComplete:[0,1,1,""],RunInfo:[0,1,1,""],RunParameters:[0,1,1,""],SampleSheetIO:[0,1,1,""],etreeToDict:[0,3,1,""],getIlluminaName:[0,3,1,""],getInfFromSeqDesc:[0,3,1,""],getInfFromSeqID:[0,3,1,""],getLibNameFromReadsPath:[0,3,1,""],getRunFolderInfo:[0,3,1,""],platformFromInstrumentSerialNumber:[0,3,1,""]},"anacore.illumina.ADSSampleSheetIO":{filterPanels:[0,2,1,""],findSplFiles:[0,2,1,""],setSplFiles:[0,2,1,""]},"anacore.illumina.SampleSheetIO":{SECTIONS:[0,4,1,""]},"anacore.maf":{MAFIO:[0,1,1,""],getName:[0,3,1,""]},"anacore.maf.MAFIO":{isRecordLine:[0,2,1,""],write:[0,2,1,""]},"anacore.matrix":{DistanceMatrixIO:[0,1,1,""]},"anacore.msi":{LocusClassifier:[0,1,1,""],LocusRes:[0,1,1,""],LocusResDistrib:[0,1,1,""],LocusResPairsCombi:[0,1,1,""],MSILocus:[0,1,1,""],MSIReport:[0,1,1,""],MSISample:[0,1,1,""],MSISplRes:[0,1,1,""],Status:[0,1,1,""],getIncompleteModels:[0,3,1,""],getNbSupporting:[0,3,1,""],toDict:[0,3,1,""]},"anacore.msi.LocusClassifier":{fit:[0,2,1,""],predict:[0,2,1,""],predict_proba:[0,2,1,""],set_status:[0,2,1,""]},"anacore.msi.LocusRes":{fromDict:[0,2,1,""]},"anacore.msi.LocusResDistrib":{fromDict:[0,2,1,""],getCount:[0,2,1,""],getDenseCount:[0,2,1,""],getDensePrct:[0,2,1,""],getMaxLength:[0,2,1,""],getMinLength:[0,2,1,""]},"anacore.msi.LocusResPairsCombi":{fromDict:[0,2,1,""],getNbFrag:[0,2,1,""]},"anacore.msi.MSILocus":{delResult:[0,2,1,""],fromDict:[0,2,1,""]},"anacore.msi.MSIReport":{parse:[0,2,1,""],write:[0,2,1,""]},"anacore.msi.MSISample":{addLocus:[0,2,1,""],delLoci:[0,2,1,""],delLocus:[0,2,1,""],fromDict:[0,2,1,""],getLociMethods:[0,2,1,""],getNbDetermined:[0,2,1,""],getNbLoci:[0,2,1,""],getNbProcessed:[0,2,1,""],getNbStable:[0,2,1,""],getNbUndetermined:[0,2,1,""],getNbUnstable:[0,2,1,""],setScore:[0,2,1,""],setStatusByInstabilityCount:[0,2,1,""],setStatusByInstabilityRatio:[0,2,1,""],setStatusByMajority:[0,2,1,""]},"anacore.msi.MSISplRes":{fromDict:[0,2,1,""]},"anacore.msi.Status":{authorizedValues:[0,2,1,""],none:[0,4,1,""],stable:[0,4,1,""],undetermined:[0,4,1,""],unstable:[0,4,1,""]},"anacore.msiannot":{MSIAnnot:[0,1,1,""],addLociResToSpl:[0,3,1,""],addLociResult:[0,3,1,""],getCastedValue:[0,3,1,""],getLocusAnnotDict:[0,3,1,""]},"anacore.msiannot.MSIAnnot":{recordToLine:[0,2,1,""]},"anacore.msings":{MSINGSAnalysis:[0,1,1,""],MSINGSReport:[0,1,1,""]},"anacore.msings.MSINGSAnalysis":{recordToLine:[0,2,1,""]},"anacore.msings.MSINGSReport":{parse:[0,2,1,""]},"anacore.node":{Node:[0,1,1,""]},"anacore.node.Node":{addChild:[0,2,1,""],fromClusterNode:[0,2,1,""],fromDict:[0,2,1,""],getAncestors:[0,2,1,""],getChildByName:[0,2,1,""],getDepth:[0,2,1,""],getDescendants:[0,2,1,""],getLeaves:[0,2,1,""],hasChild:[0,2,1,""],toDict:[0,2,1,""],toExtendedNewick:[0,2,1,""],toNewick:[0,2,1,""]},"anacore.picardIO":{PicardReader:[0,1,1,""],castCol:[0,3,1,""],getColType:[0,3,1,""],getValType:[0,3,1,""]},"anacore.region":{Region:[0,1,1,""],RegionList:[0,1,1,""],RegionTree:[0,1,1,""],consolidated:[0,3,1,""],iterOverlapped:[0,3,1,""],iterOverlappedByRegion:[0,3,1,""],mergedRegion:[0,3,1,""],splittedByRef:[0,3,1,""]},"anacore.region.Region":{contains:[0,2,1,""],getCoordinatesStr:[0,2,1,""],getMinDist:[0,2,1,""],getPosOnRef:[0,2,1,""],getPosOnRegion:[0,2,1,""],hasOverlap:[0,2,1,""],hasStrandedOverlap:[0,2,1,""],length:[0,2,1,""],setReference:[0,2,1,""],strandedContains:[0,2,1,""]},"anacore.region.RegionList":{getContainers:[0,2,1,""],getNearests:[0,2,1,""],getOverlapped:[0,2,1,""]},"anacore.region.RegionTree":{addChild:[0,2,1,""],sortChildren:[0,2,1,""]},"anacore.sequenceIO":{Faidx:[0,1,1,""],FaidxRecord:[0,1,1,""],FastaIO:[0,1,1,""],FastqIO:[0,1,1,""],IdxFastaIO:[0,1,1,""],Sequence:[0,1,1,""],SequenceFileReader:[0,1,1,""]},"anacore.sequenceIO.Faidx":{readById:[0,2,1,""]},"anacore.sequenceIO.FastaIO":{close:[0,2,1,""],isValid:[0,2,1,""],nbSeq:[0,2,1,""],nbSeqAndNt:[0,2,1,""],nextSeq:[0,2,1,""],seqToFastaLine:[0,2,1,""],write:[0,2,1,""]},"anacore.sequenceIO.FastqIO":{close:[0,2,1,""],isValid:[0,2,1,""],nbSeq:[0,2,1,""],nbSeqAndNt:[0,2,1,""],nextSeq:[0,2,1,""],qualOffset:[0,2,1,""],seqToFastqLine:[0,2,1,""],write:[0,2,1,""]},"anacore.sequenceIO.IdxFastaIO":{get:[0,2,1,""],getSub:[0,2,1,""]},"anacore.sequenceIO.Sequence":{dnaRevCom:[0,2,1,""],dna_complement:[0,4,1,""],rnaRevCom:[0,2,1,""],rna_complement:[0,4,1,""]},"anacore.sequenceIO.SequenceFileReader":{factory:[0,2,1,""]},"anacore.sv":{HashedSVIO:[0,1,1,""],SVIO:[0,1,1,""]},"anacore.sv.HashedSVIO":{recordToLine:[0,2,1,""]},"anacore.sv.SVIO":{isRecordLine:[0,2,1,""],isValid:[0,2,1,""],recordToLine:[0,2,1,""],write:[0,2,1,""]},"anacore.tophatFusion":{TopHatFusionIO:[0,1,1,""]},"anacore.vcf":{HeaderAttr:[0,1,1,""],HeaderDescAttr:[0,1,1,""],HeaderFilterAttr:[0,1,1,""],HeaderFormatAttr:[0,1,1,""],HeaderInfoAttr:[0,1,1,""],HeaderTypedAttr:[0,1,1,""],VCFIO:[0,1,1,""],VCFRecord:[0,1,1,""],decodeInfoValue:[0,3,1,""],encodeInfoValue:[0,3,1,""],getAlleleRecord:[0,3,1,""],getFreqMatrix:[0,3,1,""],getHeaderAttr:[0,3,1,""]},"anacore.vcf.HeaderAttr":{items:[0,2,1,""],keys:[0,2,1,""]},"anacore.vcf.VCFIO":{copyHeader:[0,2,1,""],isRecordLine:[0,2,1,""],recToVCFLine:[0,2,1,""],write:[0,2,1,""],writeHeader:[0,2,1,""]},"anacore.vcf.VCFRecord":{containsIndel:[0,2,1,""],fastStandardize:[0,2,1,""],getAD:[0,2,1,""],getAF:[0,2,1,""],getAFBySample:[0,2,1,""],getAltAD:[0,2,1,""],getAltAF:[0,2,1,""],getDP:[0,2,1,""],getEmptyAlleleMarker:[0,2,1,""],getMostDownstream:[0,2,1,""],getMostUpstream:[0,2,1,""],getName:[0,2,1,""],getPopAltAD:[0,2,1,""],getPopAltAF:[0,2,1,""],getPopDP:[0,2,1,""],getPopRefAD:[0,2,1,""],getPopRefAF:[0,2,1,""],isDeletion:[0,2,1,""],isIndel:[0,2,1,""],isInsertion:[0,2,1,""],normalizeSingleAllele:[0,2,1,""],refEnd:[0,2,1,""],refStart:[0,2,1,""],type:[0,2,1,""]},anacore:{STARLog:[0,0,0,"-"],abstractFile:[0,0,0,"-"],annotVcf:[0,0,0,"-"],bed:[0,0,0,"-"],filters:[0,0,0,"-"],fusion:[0,0,0,"-"],genomicRegion:[0,0,0,"-"],gff:[0,0,0,"-"],gtf:[0,0,0,"-"],hgvs:[0,0,0,"-"],illumina:[0,0,0,"-"],maf:[0,0,0,"-"],matrix:[0,0,0,"-"],msi:[0,0,0,"-"],msiannot:[0,0,0,"-"],msings:[0,0,0,"-"],node:[0,0,0,"-"],picardIO:[0,0,0,"-"],region:[0,0,0,"-"],sequenceIO:[0,0,0,"-"],sv:[0,0,0,"-"],tophatFusion:[0,0,0,"-"],vcf:[0,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","function","Python function"],"4":["py","attribute","Python attribute"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:function","4":"py:attribute"},terms:{"20ng":0,"25245351c":0,"274g":0,"3end_fusion_partn":0,"5end_fusion_partn":0,"boolean":0,"byte":0,"char":0,"class":0,"default":0,"float":0,"function":0,"int":0,"new":0,"return":0,"short":0,"static":0,"true":0,CDS:0,For:0,IDs:0,NOT:0,One:0,POS:0,The:0,These:0,Use:0,With:0,_getter:0,_parent_dist_from_root:0,_root_dist_from_leav:0,aaa:0,aaat:0,aacgt:0,aagc:0,aagt:0,aalength:0,aatc:0,absolut:0,abstractfil:1,accept_miss:0,acces:0,access:0,accord:0,account:0,acid:0,action:0,add:0,addchild:0,added:0,adding:0,addit:0,additionn:0,addlocirestospl:0,addlociresult:0,addlocu:0,addprotein:0,addtoannot:0,adssamplesheetio:0,affect:0,after:0,agc:0,age:0,age_filt:0,aggreg:0,agt:0,algorithm:0,all:0,allel:0,allow:0,alreadi:0,also:0,alt:0,altallel:0,alter:0,altern:0,alwai:0,amino:0,among:0,amplicond:0,analys:0,analysi:0,analyz:0,anapath:1,anatomo:1,ancestor:0,ani:0,anlysi:0,ann:0,ann_titl:0,annot:0,annot_field:0,annotation_field:0,annott:0,annotvcf:1,annotvcfio:0,apparit:0,appli:0,area:0,arg:0,argument:0,arrai:0,arriba:0,arribaio:0,ascend:0,atac:0,atg:0,attribut:0,author:0,authorizedvalu:0,avail:0,averag:0,base:0,basenam:0,batch:0,becom:0,bed:1,bedio:0,bedrecord:0,bedrecordtobedlin:0,been:0,befor:0,begin:0,between:0,betwen:0,bit:0,bitbucket:0,blob:0,block:0,blockcount:0,blocksiz:0,blockstart:0,bnd:0,bool:0,branch:0,breakend:0,breakendvcfio:0,breakpoint1:0,breakpoint2:0,build:0,calcul:0,call:0,callabl:0,caller:0,can:0,cancer:[0,1],cannot:0,canon:0,care:0,cast:0,cast_fct:0,castcol:0,ccb:0,cds:0,certain:0,cgt:0,chang:0,charact:0,check:0,check_ref:0,child:0,children:0,children_filt:0,children_nod:0,choos:0,chr1:0,chr3:0,chr:0,chr_po:0,chrom:0,chromosom:0,chromstart:0,cipo:0,classifi:0,clean:0,clf:0,close:0,closest_genomic_breakpoint1:0,closest_genomic_breakpoint2:0,cluster:0,clusternod:0,code:0,codon:0,codon_po:0,col_typ:0,colon:0,column:0,com:0,combin:0,come:0,comind:0,command:0,comment:0,compat:0,complement:0,complet:0,completedjobinfo:0,complex:0,compos:0,composit:0,comput:0,concern:0,confid:0,configur:0,confus:0,consensu:0,consid:0,consist:0,consolid:0,constrain:0,contain:[0,1],containsindel:0,content:0,context:0,contig:0,contigu:0,convert:0,coord:0,coordin:0,copi:0,copyhead:0,core:1,correspond:0,corrspond:0,count:0,count_msi_sampl:0,counts_of_common_mapping_read:0,coupl:0,coverage1:0,coverage2:0,creat:0,csq:0,csv:0,curr:0,curr_r:0,current:0,cytopathologi:1,data:0,data_by_locu:0,data_by_r:0,data_method_nam:0,databas:0,date_format:0,dbsnp:0,declar:0,decod:0,decodedvalu:0,decodeinfovalu:0,dedic:0,defin:0,delet:0,delloci:0,dellocu:0,delresult:0,depth:0,descend:0,describ:0,descript:0,description_data:0,descriptor:0,design:0,determin:0,dict:0,dictionari:0,differ:0,direction1:0,direction2:0,directli:0,directori:0,dirty_valu:0,discordant_m:0,displai:0,dist:0,distanc:0,distance_tag:0,distancematrixio:0,distinguish:0,distribut:0,dna:0,dna_compl:0,dnarevcom:0,doc:0,document:0,doe:0,doubl:0,downstream:0,drawn:0,dump:0,each:0,earlier:0,eaxmpl:0,ecpect:0,edu:0,effect:0,either:0,element:0,empti:0,empty_mark:0,encod:0,encoded_valu:0,encodedvalu:0,encodeinfovalu:0,end:0,end_pattern:0,ensembl:0,ensur:0,equal:0,escap:0,establish:0,estim:0,etc:0,etreetodict:0,eval:0,eval_region:0,evalu:0,evid:0,exampl:0,except:0,exclud:0,execut:0,exist:0,exon:0,exon_1_id:0,exon_2_id:0,expect:0,experimet:0,extend:0,extract:0,factori:0,fai:0,fai_path:0,faidx:0,faidxrecord:0,faiitem:0,fals:0,faq:0,faqformat:0,fasta:[0,1],fasta_format:0,fastaio:0,fastq:0,fastq_format:0,fastqio:0,faststandard:0,fcann:0,fct:0,featur:0,feature_handl:0,few:0,ffpm:0,fh_vcf:0,field:0,file:0,file_format:0,filenam:0,filepath:0,filter:1,filter_desc:0,filter_nam:0,filterpanel:0,filterscombin:0,filtersfromdict:0,findsplfil:0,finish:0,first:0,first_coord:0,first_time_to_upd:0,fisrt:0,fit:0,flag:0,flowcel:0,folder:0,follow:0,form:0,format1:0,format:0,formula:0,forward:0,found:0,fragment:0,frame:0,free:0,frequenc:0,from:0,fromclusternod:0,fromdict:0,fromstr:0,funcion:0,fusion:1,fusion_descript:0,fusion_finding_method:0,fusion_manu:0,fusion_point_for_gene_1:0,fusion_point_for_gene_2:0,fusion_sequ:0,fusion_transcript:0,fusioncatch:0,fusioncatcherio:0,fusionfileread:0,fusionnam:0,fusionvcf:0,gain:0,gdc:0,genbank:0,gene1:0,gene2:0,gene:0,gene_1_id:0,gene_1_symbol:0,gene_2_id:0,gene_2_symbol:0,gener:0,genescan:0,genom:0,genomicregion:1,get:0,getad:0,getaf:0,getafbysampl:0,getallelerecord:0,getaltad:0,getaltaf:0,getaltfromcoord:0,getancestor:0,getarea:0,getareasbychr:0,getbndinterv:0,getbyaccess:0,getcastedvalu:0,getcdsfromtranscript:0,getchildbynam:0,getcoltyp:0,getcontain:0,getcoorddictfromcoordstr:0,getcoordinatesstr:0,getcoordstr:0,getcount:0,getdensecount:0,getdenseprct:0,getdepth:0,getdescend:0,getdp:0,getemptyallelemark:0,getexectim:0,getfreqmatrix:0,getheaderattr:0,getidbynam:0,getilluminanam:0,getincompletemodel:0,getinffromseqdesc:0,getinffromseqid:0,getleav:0,getlibnamefromreadspath:0,getlocimethod:0,getlocusannotdict:0,getmaxlength:0,getmaxnbcol:0,getmindist:0,getminlength:0,getmostdownstream:0,getmostupstream:0,getnam:0,getnbdetermin:0,getnbfrag:0,getnbloci:0,getnbprocess:0,getnbstabl:0,getnbsupport:0,getnbundetermin:0,getnbunst:0,getnearest:0,getntposfromrefpo:0,getntposfromregionpo:0,getoverlap:0,getparsedrespons:0,getpopaltad:0,getpopaltaf:0,getpopdp:0,getpoprefad:0,getpoprefaf:0,getposonref:0,getposonregion:0,getprotbytr:0,getrecordvalu:0,getrequesturl:0,getrespons:0,getrunfolderinfo:0,getsortedareasbychr:0,getstrand:0,getsub:0,getsubfromrefpo:0,getsubfromregionpo:0,getter:0,getvaltyp:0,gff3:0,gff3io:0,gff3record:0,gff:1,github:0,given:0,gov:0,grandchildren:0,group:0,group_filt:0,gtf:1,gtf_path:0,gtfio:0,gzipe:0,handl:0,has:0,has_titl:0,haschild:0,hashedsvio:0,hasoverlap:0,hasstrandedoverlap:0,have:0,header:0,header_lin:0,headerattr:0,headerdescattr:0,headerfilterattr:0,headerformatattr:0,headerinfoattr:0,headertypedattr:0,hgv:1,hgvs_str:0,hierarch:0,hierarchi:0,himself:0,his:0,hiseq:0,hour:0,html:0,htslib:0,http:0,human:0,id_by_nam:0,id_to_nam:0,identifi:0,idx_alt:0,idxfastaio:0,ieee:0,ill:0,illumina:1,in_annot:0,in_b:0,in_path:0,in_report:0,includ:0,indel:0,index:[0,1],indic:0,info:0,inform:0,inherit:0,initio:0,input:0,insert:0,inspect:0,instability_threshold:0,instabl:0,instanc:0,institut:1,instrument:0,instrument_id:0,integ:0,intend:0,intern:0,introduc:0,intron:0,is_first:0,isdelet:0,isempti:0,isgzip:0,isindel:0,isinsert:0,ispredict:0,isrecordlin:0,isvalid:0,item:0,itemrgb:0,iter:0,iteroverlap:0,iteroverlappedbyregion:0,its:0,jan:0,jhu:0,job:0,json:0,junctionread:0,junctionreadcount:0,kei:0,kept:0,keyword:0,kit_i001:0,kit_v001:0,knownsnpid:0,kwarg:0,laboratoir:1,landmark:0,largeanchorsupport:0,last:0,latter:0,least:0,leav:0,leftbreakdinuc:0,leftbreakentropi:0,leftbreakpoint:0,leftgen:0,legend:0,legend_data:0,len:0,length:0,less:0,level:0,librari:[0,1],light:0,like:0,limit:0,line:0,line_bas:0,line_width:0,link:0,list:0,lite:0,load:0,loadindex:0,loadmodel:0,locat:0,loci:0,locu:0,locus_id:0,locus_nam:0,locus_weight_is_scor:0,locusclassifi:0,locusr:0,locusresdistrib:0,locusrespairscombi:0,log:0,longest_anchor_found:0,lower:0,maf:1,maf_format:0,mafio:0,manag:0,mani:0,manifest:0,manipul:0,manual:0,map:0,marker:0,master:0,matrix:1,maximum:0,memori:0,merg:0,merge_contigu:0,merge_traceback:0,merged_region:0,mergedregion:0,metadata:0,method:0,method_nam:0,microsatellit:0,million:0,min_support_model:0,min_voting_loci:0,minim:0,minimum:0,mismatch:0,miss:0,missing_replac:0,mode:0,model:0,model_method_nam:0,modulo:0,most:0,move:0,movement:0,msi:1,msi_locu:0,msi_object:0,msi_sampl:0,msi_spl:0,msiannot:1,msilocu:0,msing:1,msingsanalysi:0,msingsreport:0,msireport:0,msisampl:0,msisplr:0,mss:0,multi:0,multipl:0,must:0,mutalyz:0,mutalyzer_url:0,mutalyzerbatch:0,mutat:0,name:0,namecheck:0,nan:0,nb_loci_undetermin:0,nbseq:0,nbseqandnt:0,nc_000012:0,ncbi:0,ndaniel:0,nearest:0,necessari:0,newick:[0,1],newlin:0,next:0,nextseq:0,ng_012337:0,nm_000222:0,node:1,non:0,none:0,normal:0,normalizesingleallel:0,note:0,np_000213:0,np_001245144:0,nr21:0,nucleotid:0,number:0,numpi:0,object:0,offset:0,old:0,oncopol:1,one:0,onli:0,ontolog:0,open:0,oper:0,operon:0,opposit:0,option:0,order:0,org:0,other:0,other_nod:0,out_path:0,output:0,overlap:0,packag:1,pad:0,pair:0,panel:0,param:0,paramet:0,parent:0,parent_nod:0,pars:0,parser:0,part:0,pass:0,path:0,patient:0,pattern:0,peptide_sequ:0,per:0,percent:0,percentag:0,pfilter:0,pformat:0,phase:0,picard:0,picardio:1,picardread:0,piec:0,pip:1,placebo:0,platform:0,platformfrominstrumentserialnumb:0,popul:0,pos:0,posit:0,positionconvert:0,possibl:0,pre:0,predict:0,predict_proba:0,predicted_effect:0,previous:0,print:0,printabl:0,probabl:0,procedur:0,process:0,produc:0,progress:0,properti:0,prot:0,protein:0,protein_po:0,proteindescript:0,provid:0,proxi:0,proxy_url:0,qual:0,qual_offset:0,qualifi:0,qualiti:0,qualoffset:0,queri:0,queries_by_chr:0,quot:0,rais:0,rate:0,ratio:0,reach:0,read:0,read_identifi:0,readabl:0,readbyid:0,reader:0,reading_fram:0,readthrough:0,recfromregion:0,recommend:0,record:0,recordtolin:0,rectovcflin:0,recurs:0,reduc:0,ref:0,ref_po:0,ref_seq:0,refallel:0,refend:0,refer:0,refseq:0,refstart:0,regard:0,region:1,region_list:0,region_po:0,region_record:0,regionlist:0,regiontre:0,regul:0,reject:0,rel:0,relat:0,remov:0,replac:0,report:0,repres:0,represent:0,represet:0,request:0,requir:0,res:0,res_cl:0,respons:0,restrict_to:0,result:0,result_id:0,retri:0,retrict:0,retriev:0,retrun:0,retun:0,revers:0,rfc:0,rgb:0,rightbreakdinuc:0,rightbreakentropi:0,rightbreakpoint:0,rightgen:0,rna:0,rna_compl:0,rnarevcom:0,root:0,row:0,rtacomplet:0,rule:0,run:0,run_fold:0,runinfo:0,runmutalyz:0,runmutalyzerdescript:0,runmutalyzerlegend:0,runparamet:0,same:0,sampl:0,sample_nam:0,samplesheet:0,samplesheetio:0,sanger:0,scipi:0,score:0,sdhd_v001:0,search:0,second:0,second_coord:0,section:0,see:0,select:0,select_fct:0,selected_manifest:0,self:0,semant:0,semi:0,separ:0,seq_desc:0,seq_handl:0,seq_id:0,seq_path:0,seqnam:0,seqtofastalin:0,seqtofastqlin:0,sequenc:0,sequence_record:0,sequencefileread:0,sequenceio:1,serial:0,server:0,set:0,set_statu:0,setaggreg:0,setfct:0,setgett:0,setoper:0,setprotein:0,setrefer:0,setscor:0,setsplfil:0,setstatusbyinstabilitycount:0,setstatusbyinstabilityratio:0,setstatusbymajor:0,settranscript:0,setvcfhead:0,sever:0,shard:0,should:0,shtml:0,siblings_idx:0,sign:0,similar:0,simpli:0,singl:0,site1:0,site2:0,size:0,sjdb:0,sklearn:0,snp:0,snpconvert:0,sofa:0,softwar:0,solexa:0,sort:0,sortchildren:0,sourc:0,space:0,spanning_pair:0,spanning_unique_read:0,spanningfrag:0,spanningfragcount:0,special:0,specif:0,specifi:0,speed:0,spl_name:0,splice:0,splicetyp:0,split_reads1:0,split_reads2:0,splittedbyref:0,stabil:0,stabl:0,standard:0,star:0,starfusionio:0,starlog:1,start:0,state:0,statu:0,std:0,stop:0,store:0,str:0,strand1:0,strand2:0,strand:0,strandedcontain:0,string:0,strongli:0,structur:0,sub:0,subclass:0,subject:0,subject_by_chr:0,submiss:0,submit:0,substitut:0,success:0,suhrig:0,sum:0,superior:0,support:0,svio:0,syncrequest:0,synopsi:0,syntax:0,syntaxcheck:0,system:0,tab:0,tag:0,take:0,taken:0,target:0,term:0,test:0,test_dataset:0,text:0,than:0,thei:0,thi:0,thickend:0,thickli:0,thickstart:0,third:0,three:0,threshold:0,throughout:0,time:0,time_to_upd:0,titl:0,title_start:0,todict:0,todo:0,toextendednewick:0,togff:0,tonewick:0,too:0,tool:0,tophat:0,tophatfus:1,tophatfusionio:0,total:0,toulous:1,trace:0,track:0,train:0,train_dataset:0,transcript:0,transcript_po:0,transcriptdescript:0,treatment:0,treatment_filt:0,tree:0,tsv:0,two:0,txt:0,type:0,type_nam:0,typic:0,ucsc:0,undetermin:0,undetermined_weight:0,unecessari:0,unescap:0,uniqu:0,universitair:1,unlink:0,unmap:0,unstabl:0,unsuffici:0,updat:0,updatestatu:0,upstream:0,url:0,use:0,use_cach:0,used:0,using:0,utr:0,uwlabm:0,val:0,valid:0,valu:0,variant:0,variat:0,vcf:1,vcf_io:0,vcf_path:0,vcfio:0,vcfrecord:0,vep:0,vepvcfio:0,version:0,view:0,vote:0,wait:0,warn:0,weight:0,were:0,wheight:0,when:0,where:0,which:0,wiki:0,wikipedia:0,within:0,word:0,wrap:0,write:0,write_nb_col:0,writehead:0,www:0,xml:0,xor:0,you:0,young_with_treat:0},titles:["AnaCore package","AnaCore - managing standard file formats and objects from NGS"],titleterms:{NGS:1,abstractfil:0,anacor:[0,1],annotvcf:0,bed:0,content:1,copyright:1,descript:1,file:1,filter:0,format:1,from:1,fusion:0,genomicregion:0,gff:0,gtf:0,hgv:0,illumina:0,indic:1,instal:1,maf:0,manag:1,matrix:0,msi:0,msiannot:0,msing:0,node:0,object:1,packag:0,picardio:0,region:0,sequenceio:0,standard:1,starlog:0,tabl:1,tophatfus:0,vcf:0}})