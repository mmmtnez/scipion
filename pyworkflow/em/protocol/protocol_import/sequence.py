# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists, basename, abspath
import os
import pyworkflow.protocol.params as params
from base import ProtImportFiles
from pyworkflow.em import Sequence
from pyworkflow.em.handler_sequence import SEQ_TYPE_AMINOACIDS,\
    SEQ_TYPE_NUCLEOTIDES, IUPAC_PROTEIN_ALPHABET, SEQ_TYPE, \
    EXTENDED_PROTEIN_ALPHABET, IUPAC_NUCLEOTIDE_ALPHABET, \
    EXTENDED_DNA_ALPHABET, SequenceHandler
from pyworkflow.em. handler_atom_struct import AtomicStructHandler
from tkMessageBox import showerror

def errorWindow(tkParent, msg):
    try:
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class ProtImportSequence(ProtImportFiles):
    """ Protocol to import an aminoacid/nucleotide sequence file to the
    project"""
    _label = 'import sequence'
    SEQUENCEFILENAME = '_sequence.fasta'
    IMPORT_FROM_PLAIN_TEXT = 0
    IMPORT_FROM_STRUCTURE = 1
    IMPORT_FROM_FILES = 2
    IMPORT_FROM_UNIPROT = 3
    IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT = 0
    IMPORT_FROM_NUCLEOTIDE_FILES = 1
    IMPORT_FROM_GENEBANK = 2
    IMPORT_STRUCTURE_FROM_ID = 0
    IMPORT_STRUCTURE_FROM_FILES = 1

    url = "http://www.uniprot.org/uniprot/"

    def __init__(self, **args):
        ProtImportFiles.__init__(self, **args)

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputSequenceID', params.StringParam,
                      label="Sequence ID", allowsNull=True, important=True,
                      help="Write a sequence ID. Otherwise, if the "
                           "sequence derives from GeneBank/UniProt/PDB "
                           "databases, the respective database ID will be "
                           "selected as starting sequence ID; examples: if "
                           "you select GeneBank accession AJ520101, SCIPION "
                           "will assign AJ520101 as sequence ID; if "
                           "you select UniProt accession P12345, SCIPION will "
                           "assign P12345 as sequence ID; if you "
                           "select atomic structure 3lqd.cif, chain B, "
                           "SCIPION will assign 3lqd_B as sequence ID. In "
                           "the rest of cases, the Sequence name "
                           "will be selected as starting Sequence ID.")
        form.addParam('inputSequenceName', params.StringParam, important=True,
                      label="Sequence name", allowsNull=False,
                      help="Write a sequence name.")
        form.addParam('inputSequenceDescription', params.StringParam,
                      label="Sequence description", important=True,
                      allowsNull=True,
                      help="Write a description for your sequence. Otherwise, "
                           "if the "
                           "sequence derives from GeneBank/UniProt/PDB "
                           "databases, the respective database description "
                           "will be "
                           "selected as starting sequence description. In "
                           "the rest of cases, no sequence description will "
                           "be added.")
        form.addParam('inputSequence', params.EnumParam,
                      pointerClass='Sequence',
                      choices=SEQ_TYPE,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Import sequence of ",
                      default=SEQ_TYPE_AMINOACIDS,
                      help='Select the type of sequence to import.')
        form.addParam('inputProteinSequence', params.EnumParam,
                      choices=['plain text', 'atomic structure', 'file',
                               'UniProt ID'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d' % SEQ_TYPE_AMINOACIDS,
                      label="From ",
                      default=self.IMPORT_FROM_PLAIN_TEXT,
                      help='Select one of the four options: write the '
                           'aminoacid sequence or import it '
                           'from a previously loaded atomic structure, a local '
                           'file or an online server.')
        form.addParam('proteinIUPACalphabet', params.EnumParam,
                      choices=IUPAC_PROTEIN_ALPHABET,
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_PLAIN_TEXT),
                      label="IUPAC Protein alphabet: ",
                      default=EXTENDED_PROTEIN_ALPHABET,
                      help='Select the type of protein '
                           'alphabet:\n\nProtein alphabet: IUPAC protein '
                           'alphabet of the 20 standard amino acids; uppercase'
                           ' and single letter: *ACDEFGHIKLMNPQRSTVWY*.\n\n'
                           'Extended Protein alphabet: Extended uppercase '
                           'IUPAC '
                           'protein single letter alphabet including X etc.\n'
                           'In addition to the standard 20 single letter '
                           'protein codes, this includes:\n*B = Asx*; '
                           'Aspartic acid (R) or Asparagine (N)\n*X = Xxx*"; '
                           'Unknown or other amino acid\n*Z = Glx*; Glutamic '
                           'acid (E) or Glutamine (Q)\n*J = Xle*; Leucine ('
                           'L) or Isoleucine (I), used in mass-spec (NMR)\n'
                           '*U = Sec*; Selenocysteine\n*O = Pyl*; '
                           'Pyrrolysine\nThis alphabet is not intended to be '
                           'used with X for Selenocysteine (an ad-hoc standard'
                           ' prior to the IUPAC adoption of U instead).\n')
        form.addParam('inputRawSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_PLAIN_TEXT),
                      label="Write your sequence here:", important=True,
                      help="Write the aminoacid raw sequence.\n")
        form.addParam('inputStructureSequence', params.EnumParam,
                      choices=['id', 'file'],
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_STRUCTURE),
                      label="Atomic structure from",
                      default=self.IMPORT_STRUCTURE_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import structure data from online server or local '
                           'file',
                      pointerClass='PdbFile',
                      allowsNull=True)
        form.addParam('pdbId', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d and '
                                'inputStructureSequence == %d'
                                % (SEQ_TYPE_AMINOACIDS,
                                   self.IMPORT_FROM_STRUCTURE,
                                   self.IMPORT_STRUCTURE_FROM_ID),
                      label="Atomic structure ID ", allowsNull=True,
                      help='Type a structure ID (four alphanumeric '
                           'characters).')
        form.addParam('pdbFile', params.PathParam, label="File path",
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d and '
                                'inputStructureSequence == %d'
                                % (SEQ_TYPE_AMINOACIDS,
                                   self.IMPORT_FROM_STRUCTURE,
                                   self.IMPORT_STRUCTURE_FROM_FILES),
                      allowsNull=True,
                      help='Specify a path to desired atomic structure.')
        form.addParam('inputStructureChain', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                self.IMPORT_FROM_STRUCTURE),
                      label="Chain ", allowsNull=True,
                      help="Select a particular chain of the atomic "
                           "structure.")
        form.addParam('fileAminoacidSequence', params.PathParam,
                      label="File path",
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                self.IMPORT_FROM_FILES),
                      allowsNull=True,
                      help='Specify a path to desired aminoacid sequence '
                           'file.\nRaw sequence and fasta format are '
                           'allowed.\nIf your file contains more than one '
                           'sequence, only the first one will be considered.')
        form.addParam('uniProtSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_UNIPROT),
                      label="UniProt name/ID ", allowsNull=True,
                      help='Write a UniProt ID (six or ten alphanumeric '
                           'characters; examples: A2BC19, P12345, '
                           'A0A022YWF9, DGAL_ECOLI).\n You can convert other '
                           'database identifiers to UniProt accession codes '
                           'by using the "ID Mapping" tab on '
                           'https://www.uniprot.org/')
        form.addParam('inputNucleotideSequence', params.EnumParam,
                      choices=['plain text', 'file', 'GeneBank'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d' % SEQ_TYPE_NUCLEOTIDES,
                      label="From ",
                      default=self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT,
                      help='Select one of the four options: write the '
                           'nucleic acid sequence or import it '
                           'from a local file or an online server.')
        form.addParam('nucleotideIUPACalphabet', params.EnumParam,
                      choices=IUPAC_NUCLEOTIDE_ALPHABET,
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                 self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT),
                      label="IUPAC Nucleic acid alphabet: ",
                      default=EXTENDED_DNA_ALPHABET,
                      help='Select the type of nucleic acid '
                           'alphabet:\n\n Unambiguous DNA alphabet: '
                           'Uppercase IUPAC unambiguous DNA (letters *GATC* '
                           'only).\n\nExtended DNA: Extended IUPAC DNA '
                           'alphabet.\nIn addition to the standard letter '
                           'codes GATC, this includes:\n*B* = 5-bromouridine\n'
                           '*D* = 5,6-dihydrouridine\n*S* = thiouridine\n*W* '
                           '= wyosine\n\nAmbiguous DNA: Uppercase IUPAC '
                           'ambiguous DNA: *GATCRYWSMKHBVDN*.\n\nRNAAlphabet: '
                           'Generic single letter RNA alphabet.\n\n'
                           'Ambiguous RNA: Uppercase IUPAC ambiguous RNA; '
                           '*GAUCRYWSMKHBVDN*\n')
        form.addParam('rawNucleotideSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT),
                      label="Write your sequence here:", important=True,
                      help="Write the nucleotide raw sequence")
        form.addParam('fileNucleotideSequence', params.PathParam,
                      label="File path",
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                self.IMPORT_FROM_NUCLEOTIDE_FILES),
                      allowsNull=True,
                      help='Specify a path to desired nucleic acid sequence '
                           'file.\nRaw sequence and fasta format are '
                           'allowed.\nIf your file contains more than one '
                           'sequence, only the first one will be '
                           'considered.')
        form.addParam('geneBankSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                self.IMPORT_FROM_GENEBANK),
                      label="GeneBank accession ", allowsNull=True,
                      help='Write a GeneBank accession.\n')

    def _insertAllSteps(self):
        if self.inputSequenceID.get() is not None:
            self.id = self.inputSequenceID.get()
        else:
            self.id = ''
        self.name = self.inputSequenceName.get()
        if self.inputSequenceDescription.get() is not None:
            self.description = self.inputSequenceDescription.get()
        else:
            self.description = ''

        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            if self.inputProteinSequence == self.IMPORT_FROM_PLAIN_TEXT:
                rawSequence = self.inputRawSequence.get()
                self._insertFunctionStep('getRawSequenceStep', rawSequence)
            elif self.inputProteinSequence == self.IMPORT_FROM_STRUCTURE:
                chainId = self.inputStructureChain.get()
                self._insertFunctionStep('getSequenceOfChainStep', chainId)
            elif self.inputProteinSequence == self.IMPORT_FROM_UNIPROT:
                sequenceDB = self._getUniProtID()
                self._insertFunctionStep('sequenceDatabaseDownloadStep',
                                         sequenceDB)
            elif self.inputProteinSequence == self.IMPORT_FROM_FILES:
                self.sequenceFile = self.fileAminoacidSequence.get()
                sequenceFile = self.sequenceFile
                self._insertFunctionStep('fileDownloadStep', sequenceFile)
        else:
            if self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT:
                rawSequence = self.rawNucleotideSequence.get()
                self._insertFunctionStep('getRawSequenceStep', rawSequence)
            elif self.inputNucleotideSequence == self.IMPORT_FROM_GENEBANK:
                sequenceDB = self._getGeneBankID()
                self._insertFunctionStep('sequenceDatabaseDownloadStep',
                                         sequenceDB)
            elif self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_FILES:
                self.sequenceFile = self.fileNucleotideSequence.get()
                sequenceFile = self.sequenceFile
                self._insertFunctionStep('fileDownloadStep', sequenceFile)

        self._insertFunctionStep('createOutputStep')

    def getRawSequenceStep(self, rawSequence):
        isAminoacid = (self.inputSequence == SEQ_TYPE_AMINOACIDS)
        iUPACAlphabet = self._getAlphabet()
        seqHandler = SequenceHandler(rawSequence, iUPACAlphabet,
                                         isAminoacid=isAminoacid)
        if self.inputSequenceID.get() is not None:
            seqID = self.inputSequenceID.get()
        else:
            seqID = self.name

        if self.inputSequenceDescription.get() is not None:
            seqDescription = self.inputSequenceDescription.get()
        else:
            seqDescription = ''

        seqHandler.saveFile(self._getSeqFileName(), seqID,
                            sequence=None, name=self.name,
                            seqDescription=seqDescription)
        self.id = seqID
        self.description = seqDescription

    def getSequenceOfChainStep(self, chainId):
        self.structureHandler = AtomicStructHandler()
        if self.pdbId.get() is not None:
            pdbID = self.pdbId.get()
            inputStructure = "/tmp/" + pdbID + ".cif"
        else:
            inputStructure = self.pdbFile.get()

        self.structureHandler.read(inputStructure)
        parsed_structure = self.structureHandler.getStructure()
        self.structureHandler.editedChainSelection(parsed_structure)
        chainList = self.structureHandler.getChainList()
        chainListStr = []
        for i in chainList:
            chainListStr.append(str(i))

        if chainId in chainListStr:
            chainIndex = chainListStr.index(chainId)
            self._getSeqFileName()
            chainListSeq = self.structureHandler.getChainListSeq()
            _sequence = chainListSeq[chainIndex]
            iUPACAlphabet = self._getAlphabet()
            self.seqHandler = SequenceHandler(_sequence, iUPACAlphabet,
                                              isAminoacid=True)
            structureID = basename(inputStructure).split(".")[0].strip()
            # <class 'Bio.PDB.Chain.Chain'>
            _chain = chainList[chainIndex][1].split(":")[1].strip()
            _model = chainList[chainIndex][0].split(":")[1].strip()
            # Assignation of sequence ID: if the user has provided a specific
            #  ID, this will be adopted by default; otherwise, a sequence ID
            # related with the starting structure will be selected.
            if self.inputSequenceID.get() is not None:
                seqID = self.inputSequenceID.get()
            else:
                # Check the number of models/chains of the structure
                models = self.structureHandler.getModels()
                if len(models) > 1:
                    seqID = structureID + "_" + _model +  "_" + _chain
                else:
                    seqID = structureID + "_" + _chain
                chains = self.structureHandler.getChains()
                if len(chains) == 0:
                    seqID = seqID[:-1]

            if self.inputSequenceDescription.get() is not None:
                seqDescription = self.inputSequenceDescription.get()
            else:
                seqDescription = ''

            self.seqHandler.saveFile(self._getSeqFileName(), seqID,
                                     sequence=None, name=self.name,
                                     seqDescription=seqDescription)

            print "Selected chain: " + _chain + " from model: " + \
                  _model + " from structure: " + basename(inputStructure)
        else:
            print "Error: The selected chain ID is not included in the " \
                  "atomic structure.\n"
        self.id = seqID
        self.description = seqDescription

    def sequenceDatabaseDownloadStep(self, sequenceDB):
        """Download UniProt/GeneBank sequence from its respective database
        """
        sequenceDB = str(sequenceDB)
        if self.uniProtSequence.get() is not None:
            try:
                seqHandler = SequenceHandler()
            except Exception, e:
                error = "Error: " + sequenceDB + \
                        " is not a UniProt name/ID\n"
                errorWindow(None, error)
                self.protocol.setAborted()
        elif self._getGeneBankID() is not None:
            try:
                seqHandler = SequenceHandler(isAminoacid=False)
            except Exception, e:
                error = "Error: " + str(sequenceDB) + \
                        " is not a GeneBank name/ID\n"
                errorWindow(None, error)
                self.protocol.setAborted()

        record = seqHandler.downloadSeqFromDatabase(sequenceDB)

        if self.inputSequenceID.get() is not None:
            seqID = self.inputSequenceID.get()
        elif sequenceDB != '':
            seqID = sequenceDB
        else:
            seqID = self.name
        if self.inputSequenceDescription.get() is not None:
            seqDescription = self.inputSequenceDescription.get()
        elif record.description != '':
            seqDescription = record.description
        else:
            seqDescription = ''

        seqFileName = self._getSeqFileName()
        seqHandler.saveFile(seqFileName, seqID, sequence=record.seq,
                            name=self.name,
                            seqDescription=seqDescription)
        self.id = seqID
        self.description = seqDescription

    def fileDownloadStep(self, sequenceFile):
        # If sequencePath contains more than one sequence, only
        # the first one will be considered
        seqHandler = SequenceHandler()
        record = seqHandler.downloadSeqFromFile(sequenceFile,
                                                        type="fasta")
        if self.inputSequenceID.get() is not None:
            seqID = self.inputSequenceID.get()
        elif record.id != '':
            seqID = record.id
        else:
            seqID = self.name
        if self.inputSequenceDescription.get() is not None:
            seqDescription = self.inputSequenceDescription.get()
        elif record.description != '':
            seqDescription = record.description
        else:
            seqDescription = ''
        seqFileName = self._getSeqFileName()
        seqHandler.saveFile(seqFileName, seqID, sequence=record.seq,
                            name=self.name,
                            seqDescription=seqDescription)
        self.id = seqID
        self.description = seqDescription

    def createOutputStep(self):
        """ Register the output object. """
        seq = Sequence(name=self.name, filename=self._getSeqFileName(),
                       isAminoacids=(self.inputSequence ==
                                     SEQ_TYPE_AMINOACIDS),
                       id=self.id, description=self.description)
        outputs = {'outputSequence': seq}
        self._defineOutputs(**outputs)


    def _summary(self):
        summary = []
        self.name = self.inputSequenceName.get()
        #inputStructure = self._getStructureName()
        uniProtId = self._getUniProtID()
        geneBankID = self._getGeneBankID()
        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            summary.append('Sequence of aminoacids:\n')
            if self.inputProteinSequence == self.IMPORT_FROM_PLAIN_TEXT:
                summary.append("Sequence *%s* imported from plain text\n"
                               % self.name)
            elif self.inputProteinSequence == self.IMPORT_FROM_STRUCTURE:
                if self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_ID:
                    summary.append("Sequence *%s* imported from atomic "
                                   "structure *%s.cif*\n"
                                   % (self.name, self.pdbId.get()))
                elif self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_FILES:
                    summary.append("Sequence *%s* imported from file *%s*\n"
                                   % (self.name, self.pdbFile.get()))
            elif self.inputProteinSequence == self.IMPORT_FROM_UNIPROT:
                summary.append("Sequence *%s* imported from UniProt ID "
                               "*%s*\n"
                               % (self.name, uniProtId))
            elif self.inputProteinSequence == self.IMPORT_FROM_FILES:
                summary.append("Sequence *%s* imported from file name: "
                               "*%s*\n"
                               % (self.name, self.fileAminoacidSequence.get()))
        else:
            summary.append('Sequence of nucleotides:\n')
            if self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT:
                summary.append("Sequence *%s* imported from plain text\n"
                               % self.name)
            elif self.inputNucleotideSequence == self.IMPORT_FROM_GENEBANK:
                summary.append("Sequence *%s* imported from geneBank ID "
                               "*%s*\n"
                               % (self.name, geneBankID))
            elif self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_FILES:
                summary.append("Sequence *%s* imported from file name: "
                               "*%s*\n"
                               % (self.name, self.fileNucleotideSequence.get()))
        return summary

    def _validate(self):
        errors = []
        return errors

    def _getSequenceName(self):
        pass

    def _getSeqFileName(self):
        SEQUENCEFILENAME = self._getExtraPath(self.name + self.SEQUENCEFILENAME)
        return os.path.abspath(SEQUENCEFILENAME)

    def _getFastaFormat(self, seqID):
        self.fasta = ">" + seqID

    def _getUniProtID(self):
        return self.uniProtSequence.get()

    def _getGeneBankID(self):
        return self.geneBankSequence

    def _getAlphabet(self):
        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            return self.proteinIUPACalphabet.get()
        else:
            return self.nucleotideIUPACalphabet.get()

