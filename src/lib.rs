#![doc = include_str!("../README.md")]

extern crate purr;

use std::collections::HashMap;

use purr::feature::Aliphatic;
use purr::feature::Aromatic;
use purr::feature::AtomKind;
use purr::feature::BondKind;
use purr::feature::BracketSymbol;
use purr::feature::Configuration;
use purr::feature::Element;
use purr::feature::VirtualHydrogen;
use purr::feature::Rnum;
use purr::walk::Follower;

/// An error marker for sequences containing invalid amino acids.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnknownResidue;

impl std::fmt::Display for UnknownResidue {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        write!(f, "unknown residue found in sequence")
    }
}

impl std::error::Error for UnknownResidue {}

/// A generic error type for this crate.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// A cross-link is invalid.
    ///
    /// This issue can occur when a requested cross-link cannot be created
    /// from the amino acid residues at the given locations.
    ///
    /// # Example
    /// A disulfide bond cannot be created from two L-alanine residues:
    /// ```rust
    /// use proteinogenic::Error;
    /// use proteinogenic::AminoAcid::Ala;
    ///
    /// let mut prot = proteinogenic::Protein::new([Ala, Ala]);
    /// prot.cross_link(proteinogenic::CrossLink::Cystine(1, 2)).unwrap();
    ///
    /// let mut f = purr::write::Writer::new();
    /// assert!(matches!(prot.visit(&mut f), Err(Error::InvalidCrossLink(1, _, _))));
    ///
    /// ```
    InvalidCrossLink(u16, AminoAcid, CrossLink),

    /// A residue is involved in more than one cross-link.
    ///
    /// # Example
    /// A cysteine can only be involved in one disulfide bond:
    /// ```rust
    /// use proteinogenic::Error;
    /// use proteinogenic::AminoAcid::Cys;
    ///
    /// let mut prot = proteinogenic::Protein::new([Cys, Cys, Cys]);
    /// prot.cross_link(proteinogenic::CrossLink::Cystine(1, 2)).unwrap();
    /// assert_eq!(
    ///     prot.cross_link(proteinogenic::CrossLink::Cystine(1, 3)),
    ///     Err(Error::DuplicateCrossLink(1)),
    /// );
    /// ```
    DuplicateCrossLink(u16),

    /// Too many cross-links were created.
    ///
    /// This can occur when a protein contains too many cross-links, which will
    /// exhaust the number of possibilites for ring identifiers in SMILES.
    TooManyCrossLinks,
}

/// A single L-α amino-acid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AminoAcid {
    /// [L-arginine](https://en.wikipedia.org/wiki/Arginine).
    ///
    /// ![Skeletal formula of L-arginine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29952)
    Arg,

    /// [L-histidine](https://en.wikipedia.org/wiki/Histidine).
    ///
    /// ![Skeletal formula of L-histidine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29979)
    His,

    /// [L-lysine](https://en.wikipedia.org/wiki/Lysine).
    ///
    /// ![Skeletal formula of L-lysine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29967)
    Lys,

    /// [L-aspartic acid (L-aspartate)](https://en.wikipedia.org/wiki/Aspartic_acid).
    ///
    /// ![Skeletal formula of L-aspartic acid](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29958)
    Asp,

    /// [L-glutamic acid (L-glutamate)](https://en.wikipedia.org/wiki/Glutamic_acid).
    ///
    /// ![Skeletal formula of L-glutamic acid](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29972)
    Glu,

    /// [L-serine](https://en.wikipedia.org/wiki/Serine).
    ///
    /// ![Skeletal formula of L-serine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29999)
    Ser,

    /// [L-threonine](https://en.wikipedia.org/wiki/Threonine).
    ///
    /// ![Skeletal formula of L-threonine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30013)
    Thr,

    /// [L-asparagine](https://en.wikipedia.org/wiki/Asparagine).
    ///
    /// ![Skeletal formula of L-asparagine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=50347)
    Asn,

    /// [L-glutamine](https://en.wikipedia.org/wiki/Glutamine).
    ///
    /// ![Skeletal formula of L-glutamine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30011)
    Gln,

    /// [Glycine](https://en.wikipedia.org/wiki/Glycine).
    ///
    /// ![Skeletal formula of Glycine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29947).
    Gly,

    /// [L-proline](https://en.wikipedia.org/wiki/Proline).
    ///
    /// ![Skeletal formula of L-proline](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=50342)
    Pro,

    /// [L-cysteine](https://en.wikipedia.org/wiki/Cysteine).
    ///
    /// ![Skeletal formula of L-cysteine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29950)
    Cys,

    /// [L-selenocysteine](https://en.wikipedia.org/wiki/Selenocysteine).
    ///
    /// ![Skeletal formula of L-selenocysteine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30000)
    Sec,

    /// [L-alanine](https://en.wikipedia.org/wiki/Alanine).
    ///
    /// ![Skeletal formula of L-alanine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=46217)
    Ala,

    /// [L-valine](https://en.wikipedia.org/wiki/Valine).
    ///
    /// ![Skeletal formula of L-valine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30015)
    Val,

    /// [L-isoleucine](https://en.wikipedia.org/wiki/Isoleucine).
    ///
    /// ![Skeletal formula of L-isoleucine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30009)
    Ile,

    /// [L-leucine](https://en.wikipedia.org/wiki/Leucine).
    ///
    /// ![Skeletal formula of L-leucine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=30006)
    Leu,

    /// [L-methionine](https://en.wikipedia.org/wiki/Methionine).
    ///
    /// ![Skeletal formula of L-methionine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=16044)
    Met,

    /// [L-phenylalanine](https://en.wikipedia.org/wiki/Phenylalanine).
    ///
    /// ![Skeletal formula of L-phenylalanine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29997)
    Phe,

    /// [L-tyrosine](https://en.wikipedia.org/wiki/Tyrosine).
    ///
    /// ![Skeletal formula of L-tyrosine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=46858)
    Tyr,

    /// [L-tryptophan](https://en.wikipedia.org/wiki/Tryptophan).
    ///
    /// ![Skeletal formula of L-tryptophan](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=29954)
    Trp,

    /// [L-Pyrrolysine](https://fr.wikipedia.org/wiki/Pyrrolysine).
    ///
    /// ![Skeletal formula of L-pyrrolysine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=21860)
    Pyl,

    /// [2,3-didehydroalanine](https://en.wikipedia.org/wiki/Dehydroalanine).
    ///
    /// ![Skeletal formula of 2,3-didehydroalanine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=90873)
    Dha,

    /// (Z)-dehydrobutyrine.
    ///
    /// ![Skeletal formula of (Z)-dehydrobutyrine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=18820)
    Dhb,
}

impl AminoAcid {
    /// Create an `AminoAcid` variant from a 1-letter code.
    pub fn from_code1(code: char) -> Result<AminoAcid, UnknownResidue> {
        match code {
            'R' => Ok(AminoAcid::Arg),
            'H' => Ok(AminoAcid::His),
            'K' => Ok(AminoAcid::Lys),
            'D' => Ok(AminoAcid::Asp),
            'E' => Ok(AminoAcid::Glu),
            'S' => Ok(AminoAcid::Ser),
            'T' => Ok(AminoAcid::Thr),
            'N' => Ok(AminoAcid::Asn),
            'Q' => Ok(AminoAcid::Gln),
            'G' => Ok(AminoAcid::Gly),
            'P' => Ok(AminoAcid::Pro),
            'C' => Ok(AminoAcid::Cys),
            'U' => Ok(AminoAcid::Sec),
            'A' => Ok(AminoAcid::Ala),
            'V' => Ok(AminoAcid::Val),
            'I' => Ok(AminoAcid::Ile),
            'L' => Ok(AminoAcid::Leu),
            'M' => Ok(AminoAcid::Met),
            'F' => Ok(AminoAcid::Phe),
            'Y' => Ok(AminoAcid::Tyr),
            'W' => Ok(AminoAcid::Trp),
            'O' => Ok(AminoAcid::Pyl),
            _ => Err(UnknownResidue),
        }
    }

    /// Create an `AminoAcid` variant from a 3-letter code.
    pub fn from_code3(code: &str) -> Result<AminoAcid, UnknownResidue> {
        match code {
            "Arg" => Ok(AminoAcid::Arg),
            "His" => Ok(AminoAcid::His),
            "Lys" => Ok(AminoAcid::Lys),
            "Asp" => Ok(AminoAcid::Asp),
            "Glu" => Ok(AminoAcid::Glu),
            "Ser" => Ok(AminoAcid::Ser),
            "Thr" => Ok(AminoAcid::Thr),
            "Asn" => Ok(AminoAcid::Asn),
            "Gln" => Ok(AminoAcid::Gln),
            "Gly" => Ok(AminoAcid::Gly),
            "Pro" => Ok(AminoAcid::Pro),
            "Cys" => Ok(AminoAcid::Cys),
            "Sec" => Ok(AminoAcid::Sec),
            "Ala" => Ok(AminoAcid::Ala),
            "Val" => Ok(AminoAcid::Val),
            "Ile" => Ok(AminoAcid::Ile),
            "Leu" => Ok(AminoAcid::Leu),
            "Met" => Ok(AminoAcid::Met),
            "Phe" => Ok(AminoAcid::Phe),
            "Tyr" => Ok(AminoAcid::Tyr),
            "Trp" => Ok(AminoAcid::Trp),
            "Pyl" => Ok(AminoAcid::Pyl),
            "Dha" => Ok(AminoAcid::Dha),
            "Dhb" => Ok(AminoAcid::Dhb),
            _ => Err(UnknownResidue),
        }
    }
}

/// A covalent bond between several amino-acid residues.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CrossLink {
    /// [L-cystine](https://en.wikipedia.org/wiki/Cystine).
    ///
    /// ![Skeletal formula of L-cystine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=50058)
    Cystine(u16, u16),

    /// [meso-lanthionine](https://en.wikipedia.org/wiki/Lanthionine).
    ///
    /// ![Skeletal formula of meso-lanthionine](https://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=25205)
    Lan(u16, u16),

    /// β-methyllanthionine.
    MeLan(u16, u16),
}

/// A peptide cyclization mechanism.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cyclization {
    /// No cyclization, resulting in a linear peptide.
    None,

    /// Head-to-tail cyclization, resulting in an homodetic cyclic peptide.
    HeadToTail,
}

impl Default for Cyclization {
    fn default() -> Self {
        Cyclization::None
    }
}

/// A protein abstracted as a modified peptide.
#[derive(Debug, Clone, PartialEq)]
pub struct Protein<S> {
    cyclization: Cyclization,

    cross_links: HashMap<u16, (Rnum, CrossLink)>,
    cross_link_num: u16,

    sequence: S,
}

impl<S> Protein<S> {
    /// Mark whether the peptide is cyclized through a known cyclization mechanism.
    pub fn cyclization(&mut self, cyclization: Cyclization) -> &mut Self {
        self.cyclization = cyclization;
        self
    }

    /// Add a cross-link between residues of the peptide.
    pub fn cross_link(&mut self, cross_link: CrossLink) -> Result<&mut Self, Error> {
        let rnum = Rnum::try_from( self.cross_link_num ).unwrap(); // FIXME
        match cross_link {
            CrossLink::Cystine(i, j) | CrossLink::Lan(i, j) | CrossLink::MeLan(i, j) => {
                let val = (rnum, cross_link);
                if let Some(_) = self.cross_links.insert(i, val.clone()) {
                    return Err(Error::DuplicateCrossLink(i));
                }
                if let Some(_) = self.cross_links.insert(j, val) {
                    return Err(Error::DuplicateCrossLink(j));
                }
            }
        }

        self.cross_link_num += 1;
        Ok(self)
    }

    /// Perform a walk on the atoms and bonds of the amino acid.
    ///
    /// The follower must have been initialized with a head. It will finish its
    /// walk on the β carbon, without visiting the atoms part of the peptidic
    /// bond.
    ///
    /// # Note
    /// This is implemented as an associated method rather than a function to
    /// make the borrow-checker happy about us borrowing `self.sequence` mutably
    /// and `self.cross_links` immutably.
    ///
    fn visit_residue<F: Follower>(
        aa: AminoAcid,
        follower: &mut F,
        index: u16,
        cross_links: &HashMap<u16, (Rnum, CrossLink)>
    ) -> Result<(), Error> {
        const CARBON_TH2: AtomKind = AtomKind::Bracket {
            symbol: BracketSymbol::Element(Element::C),
            configuration: Some(Configuration::TH2),
            hcount: Some(VirtualHydrogen::H1),
            isotope: None,
            charge: None,
            map: None,
        };
        const CARBON_TH1: AtomKind = AtomKind::Bracket {
            symbol: BracketSymbol::Element(Element::C),
            configuration: Some(Configuration::TH1),
            hcount: Some(VirtualHydrogen::H1),
            isotope: None,
            charge: None,
            map: None,
        };

        // only L-threonine and L-cysteine can build a cross-link at the
        // moment, any other amino-acid has to be an error.
        if aa != AminoAcid::Thr && aa != AminoAcid::Cys {
            if let Some((_, cross_link)) = cross_links.get(&index) {
                return Err(Error::InvalidCrossLink(index, aa, *cross_link));
            }
        }

        // visit the alpha carbon and the residue
        match aa {
            AminoAcid::Dhb => {
                // alpha carbon
                follower.extend(BondKind::Up, AtomKind::Aliphatic(Aliphatic::C));
                // residue
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Down, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(2);
            }

            AminoAcid::Dha => {
                // alpha carbon
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                // residue
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
            }

            AminoAcid::Pyl => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, CARBON_TH1);
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, CARBON_TH1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::N));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.pop(11);
            }

            AminoAcid::Gly => {
                // alpha carbon
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
            }

            AminoAcid::Ala => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
            }

            AminoAcid::Pro => {
                // proline ring
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH1);
                follower.join(BondKind::Elided, Rnum::R1);
            }

            AminoAcid::Val => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(2);
            }

            AminoAcid::Leu => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(3);
            }

            AminoAcid::Met => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(4);
            }

            AminoAcid::Phe => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.pop(7);
            }

            AminoAcid::Tyr => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.pop(7);
            }

            AminoAcid::Cys => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                match cross_links.get(&index) {
                    // no cross-link, just add the thiol group.
                    None => {
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                        follower.pop(2);
                    }
                    // cystine, add the first sulfur, the other Cys will add the second one.
                    Some((rnum, CrossLink::Cystine(_, _))) => {
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                        follower.join(BondKind::Elided, rnum.clone());
                        follower.pop(2);
                    }
                    // lanthionine, by convention the sulfur comes from the first residue
                    Some((rnum, CrossLink::Lan(i, _))) if *i == index => {
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                        follower.join(BondKind::Elided, rnum.clone());
                        follower.pop(2);
                    }
                    Some((rnum, CrossLink::Lan(_, _))) => {
                        follower.join(BondKind::Elided, rnum.clone());
                        follower.pop(1);
                    }
                    // methyllanthionine, add the sulfur, the threonine won't add the hydroxy group
                    Some((rnum, CrossLink::MeLan(_, _))) => {
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                        follower.join(BondKind::Elided, rnum.clone());
                        follower.pop(2);
                    }
                }
            }

            AminoAcid::Ser => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(2);
            }

            AminoAcid::Sec => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(
                    BondKind::Elided,
                    AtomKind::Bracket {
                        symbol: BracketSymbol::Element(Element::Se),
                        isotope: None,
                        configuration: None,
                        hcount: None,
                        charge: None,
                        map: None,
                    },
                );
                follower.pop(2);
            }

            AminoAcid::Thr => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, CARBON_TH2);
                match cross_links.get(&index) {
                    None => {
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                        follower.pop(1);
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                        follower.pop(2);
                    }
                    Some((rnum, CrossLink::MeLan(_, _))) => {
                        follower.join(BondKind::Elided, rnum.clone());
                        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                        follower.pop(2);
                    }
                    Some((_, other)) => {
                        return Err(Error::InvalidCrossLink(index, aa, *other));
                    }
                }
            }

            AminoAcid::Asn => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.pop(3);
            }

            AminoAcid::Gln => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.pop(4);
            }

            AminoAcid::Arg => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::N));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.pop(6);
            }

            AminoAcid::Lys => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.pop(5);
            }

            AminoAcid::His => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.pop(6);
            }

            AminoAcid::Asp => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(3);
            }

            AminoAcid::Glu => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(4);
            }

            AminoAcid::Ile => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, CARBON_TH2);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(3);
            }

            AminoAcid::Trp => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R2);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, Rnum::R2);
                follower.pop(10);
            }
        }

        // visit the beta carbon and finish
        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
        Ok(())
    }
}

impl<S> Protein<S>
where
    S: IntoIterator<Item=AminoAcid>,
{
    pub fn new(sequence: S) -> Self {
        Self {
            sequence,
            cyclization: Cyclization::default(),
            cross_links: HashMap::new(),
            cross_link_num: 3, // R0 is used for cyclization, R1 and R2 in residues
        }
    }

    pub fn visit<F: Follower>(self, follower: &mut F) -> Result<(), Error> {
        // visit every amino acid one by one
        let mut aa_iter = self.sequence.into_iter().enumerate();
        if let Some((index, aa)) = aa_iter.next() {
            // N-terminus: create a the N of the primary amine.
            follower.root(AtomKind::Aliphatic(Aliphatic::N));
            if self.cyclization == Cyclization::HeadToTail {
                follower.join(BondKind::Elided, Rnum::R0);
            }

            // visit residue
            Self::visit_residue(aa, follower, index as u16 + 1, &self.cross_links)?;

            // add the carboxy group to the β carbon.
            follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
            follower.pop(1);
            // keep visiting following amino acids.
            while let Some((index, aa)) = aa_iter.next() {
                // next amino acid: create the N atom of the carboxamide and visit residue.
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
                Self::visit_residue(aa, follower, index as u16 + 1, &self.cross_links)?;
                // add the carboxy group to the β carbon.
                follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
            }

            // C-terminus: create the O atom of the carboxylic acid.
            if self.cyclization == Cyclization::HeadToTail {
                follower.join(BondKind::Elided, Rnum::R0);
            } else {
                follower.extend(BondKind::Single, AtomKind::Aliphatic(Aliphatic::O));
            }
        }

        Ok(())
    }
}

/// Perform a walk on the atoms and bonds of the protein.
pub fn visit<'aa, S, F>(sequence: S, follower: &mut F) -> Result<(), Error>
where
    S: IntoIterator<Item = AminoAcid>,
    F: Follower,
{
    Protein::new(sequence).visit(follower)
}

/// Create a SMILES string for the given amino-acid sequence.
pub fn smiles<'aa, S>(sequence: S) -> Result<String, Error>
where
    S: IntoIterator<Item = AminoAcid>,
{
    let mut writer = purr::write::Writer::new();
    visit(sequence, &mut writer)?;
    Ok(writer.write())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn empty() {
        let s = smiles([]).unwrap();
        assert_eq!(s, "");
    }

    #[test]
    fn from_code1() {
        assert_eq!(AminoAcid::from_code1('Y'), Ok(AminoAcid::Tyr));
        assert_eq!(AminoAcid::from_code1('α'), Err(UnknownResidue));
    }

    #[test]
    fn from_code3() {
        assert_eq!(AminoAcid::from_code3("Thr"), Ok(AminoAcid::Thr));
        assert_eq!(AminoAcid::from_code3("Xyz"), Err(UnknownResidue));
    }
}
