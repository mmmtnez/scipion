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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol.protocol_import.sequence import \
    ProtImportSequence
from pyworkflow.em.handler_sequence import \
    SEQ_TYPE_NUCLEOTIDES, PROTEIN_ALPHABET, \
    AMBIGOUS_RNA_ALPHABET
import os
from Bio import SeqIO

class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportSequence(TestImportBase):
    USERID = "UserID"
    NAME = 'USER_SEQ'
    DESCRIPTION = 'User description'
    NUCLEOTIDESEQ1 = 'AATGCGGTTGGGBDSW********GGCACACG'
    AMINOACIDSSEQ1 = 'LARKJLAKPABXZJUO********VAVAVALK'
    CHAIN = "['model: 0', 'chain: B', '146 residues']"
    pdbID = "3lqd"
    GENEBANKID = 'AJ520101.1'
    UNIPROTID = 'P12345'

    def testImportUserNucleotideSequence1(self):
        """
        Import a single nucleotide sequence provided by the user (nucleotide
        alphabet by default)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'rawNucleotideSequence': self.NUCLEOTIDESEQ1
               }
        prot1 = self.newProtocol(ProtImportSequence, **args)
        prot1.setObjLabel('1_import nucleotide,\nseq from user\nExtended DNA '
                           'alphabet')
        self.launchProtocol(prot1)
        sequence = prot1.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("USER_SEQ" in record.id)
        self.assertTrue("USER_SEQ" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("AATGCGGTTGGGBDSWGGCACACG" in record.seq)

    def testImportUserNucleotideSequence2(self):
        """
        Import a single nucleotide sequence provided by the user (nucleotide
        alphabet AMBIGOUS_RNA_ALPHABET)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'nucleotideIUPACalphabet': AMBIGOUS_RNA_ALPHABET,
                'rawNucleotideSequence': self.NUCLEOTIDESEQ1
               }
        prot2 = self.newProtocol(ProtImportSequence, **args)
        prot2.setObjLabel('2_import nucleotide,\nseq from user\nAmbigous RNA '
                           'alphabet')
        self.launchProtocol(prot2)
        sequence = prot2.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("USER_SEQ" in record.id)
        self.assertTrue("USER_SEQ" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("AAGCGGGGGBDSWGGCACACG" in record.seq)

    def testImportFileNucleotideSequence1(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileNucleotideSequence': self.dsModBuild.getFile(
                    'Sequences/AJ520101.fasta')
                }
        prot3 = self.newProtocol(ProtImportSequence, **args)
        prot3.setObjLabel('3_import nucleotide,\nseq from file')
        self.launchProtocol(prot3)
        sequence = prot3.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("AJ520101.1" in record.id)
        self.assertTrue("AJ520101.1" in record.name)
        self.assertTrue("Rhizobium leguminosarum bv. viciae plasmid "
                        "partial fixC gene, fixX gene, nifA gene "
                        "and nifB gene" in record.description)
        self.assertTrue("GGATCCGAGA" in record.seq[:10])

    def testImportFileNucleotideSequence2(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileNucleotideSequence': self.dsModBuild.getFile(
                    'Sequences/AJ520101.fasta')
                }
        prot4 = self.newProtocol(ProtImportSequence, **args)
        prot4.setObjLabel('4_import nucleotide,\nseq from file')
        self.launchProtocol(prot4)
        sequence = prot4.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("GGATCCGAGA" in record.seq[:10])

    def testImportFileNucleotideSequence3(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileNucleotideSequence': self.dsModBuild.getFile(
                    'Sequences/Several_sequences.fasta')
                }
        prot4_1 = self.newProtocol(ProtImportSequence, **args)
        prot4_1.setObjLabel('4_1_import nucleotide,\nseq from file')
        self.launchProtocol(prot4_1)
        sequence = prot4_1.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("Seq1" in record.id)
        self.assertTrue("Seq1" in record.name)
        self.assertTrue("This is the sequence number one"
                        in record.description)
        self.assertTrue("TGGCTAAATA" in record.seq[:10])

    def testImportFileNucleotideSequence4(self):
        """
        Import a single nucleotide sequence from a text file that contains
        several sequences
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileNucleotideSequence': self.dsModBuild.getFile(
                    'Sequences/Several_sequences.fasta')
                }
        prot4_2 = self.newProtocol(ProtImportSequence, **args)
        prot4_2.setObjLabel('4.2_import nucleotide,\nseq from file')
        self.launchProtocol(prot4_2)
        sequence = prot4_2.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("TGGCTAAATA" in record.seq[:10])

    def testImportGeneBankNucleotideSequence1(self):
        """
        Import a single nucleotide sequence from GeneBank
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_GENEBANK,
                'geneBankSequence': self.GENEBANKID
                }
        prot5 = self.newProtocol(ProtImportSequence, **args)
        prot5.setObjLabel('5_import nucleotide,\nseq from '
                           'GeneBank')
        self.launchProtocol(prot5)
        sequence = prot5.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("AJ520101.1" in record.id)
        self.assertTrue("AJ520101.1" in record.name)
        self.assertTrue("Rhizobium leguminosarum bv. viciae plasmid "
                        "partial fixC gene, fixX gene, nifA gene "
                        "and nifB gene" in record.description)
        self.assertTrue("GGATCCGAGA" in record.seq[:10])

    def testImportGeneBankNucleotideSequence2(self):
        """
        Import a single nucleotide sequence from GeneBank
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_GENEBANK,
                'geneBankSequence': self.GENEBANKID
                }
        prot6 = self.newProtocol(ProtImportSequence, **args)
        prot6.setObjLabel('6_import nucleotide,\nseq from '
                           'GeneBank')
        self.launchProtocol(prot6)
        sequence = prot6.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("GGATCCGAGA" in record.seq)

    def testImportUserAminoacidSequence1(self):
        """
        Import a single aminoacid sequence provided by the user (protein
        alphabet by default)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputRawSequence': self.AMINOACIDSSEQ1
                }
        prot7 = self.newProtocol(ProtImportSequence, **args)
        prot7.setObjLabel('7_import aminoacid seq,\n from user\nExtended '
                           'protein alphabet')
        self.launchProtocol(prot7)

        sequence = prot7.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("USER_SEQ" in record.id)
        self.assertTrue("USER_SEQ" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("LARKJLAKPA" in record.seq[:10])

    def testImportUserAminoacidSequence2(self):
        """
        Import a single aminoacid sequence provided by the user (protein
        alphabet PROTEIN_ALPHABET
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'proteinIUPACalphabet': PROTEIN_ALPHABET,
                'inputRawSequence': self.AMINOACIDSSEQ1
                }
        prot8 = self.newProtocol(ProtImportSequence, **args)
        prot8.setObjLabel('8_import aminoacid seq,\n from user\nProtein '
                           'alphabet')
        self.launchProtocol(prot8)
        sequence = prot8.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("USER_SEQ" in record.id)
        self.assertTrue("USER_SEQ" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("LARKLAKPAVAVAVALK" in record.seq)

    def testImportStructureAminoacidSequence1(self):
        """
        Import the sequence of chain B of atomic structure 3lqd.cif
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'proteinIUPACalphabet': PROTEIN_ALPHABET,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID,
                'inputStructureChain': self.CHAIN
                }
        prot9 = self.newProtocol(ProtImportSequence, **args)
        prot9.setObjLabel('9_import aminoacid seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot9)
        sequence = prot9.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("3lqd_B" in record.id)
        self.assertTrue("3lqd_B" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("VHLSGEEKSA" in record.seq[:10])

    def testImportStructureAminoacidSequence2(self):
        """
        Import the sequence of chain B of atomic structure 3lqd.cif
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/3lqd.cif'),
                'inputStructureChain': self.CHAIN
                }
        prot10 = self.newProtocol(ProtImportSequence, **args)
        prot10.setObjLabel('10_import aminoacid seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot10)
        sequence = prot10.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description" in record.description)
        self.assertTrue("VHLSGEEKSA" in record.seq[:10])

    def testImportFileAminoacidSequence1(self):
        """
        Import a single aminoacid sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileAminoacidSequence': self.dsModBuild.getFile(
                    'Sequences/COX1_human.fasta')
                }
        prot11 = self.newProtocol(ProtImportSequence, **args)
        prot11.setObjLabel('11_import aminoacid,\nseq from file')
        self.launchProtocol(prot11)
        sequence = prot11.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("YP_003024028.1" in record.id)
        self.assertTrue("YP_003024028.1" in record.name)
        self.assertTrue("cytochrome c oxidase subunit I (mitochondrion) "
                        "[Homo sapiens]"
                        in record.description)
        self.assertTrue("MFADRWLFST" in record.seq[:10])

    def testImportFileAminoacidSequence2(self):
        """
        Import a single aminoacid sequence from a text file
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileAminoacidSequence': self.dsModBuild.getFile(
                    'Sequences/COX1_human.fasta')
                }
        prot12 = self.newProtocol(ProtImportSequence, **args)
        prot12.setObjLabel('12_import aminoacid,\nseq from file')
        self.launchProtocol(prot12)
        sequence = prot12.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description"
                        in record.description)
        self.assertTrue("MFADRWLFST" in record.seq[:10])

    def testImportUniprotAminoacidSequence1(self):
        """
        Import a single aminoacid sequence from UniProt
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID
                }
        prot13 = self.newProtocol(ProtImportSequence, **args)
        prot13.setObjLabel('13_import aminoacids,\nseq from '
                           'UniProt')
        self.launchProtocol(prot13)
        sequence = prot13.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("P12345" in record.id)
        self.assertTrue("P12345" in record.name)
        self.assertTrue("Aspartate aminotransferase, mitochondrial"
                        in record.description)
        self.assertTrue("MALLHSARVL" in record.seq[:10])

    def testImportUniprotAminoacidSequence2(self):
        """
        Import a single aminoacid sequence from UniProt
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID
                    }
        prot14 = self.newProtocol(ProtImportSequence, **args)
        prot14.setObjLabel('14_import aminoacids,\nseq from '
                            'UniProt')
        self.launchProtocol(prot14)
        sequence = prot14.outputSequence
        seqFileName = sequence.getFileName()
        self.assertEqual(os.path.basename(seqFileName),
                         'USER_SEQ_sequence.fasta')
        record = SeqIO.read(seqFileName, "fasta")
        self.assertTrue("UserID" in record.id)
        self.assertTrue("UserID" in record.name)
        self.assertTrue("User description"
                        in record.description)
        self.assertTrue("MALLHSARVL" in record.seq[:10])


