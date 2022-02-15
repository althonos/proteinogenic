extern crate purr;

use purr::feature::Aliphatic;
use purr::feature::Aromatic;
use purr::feature::AtomKind;
use purr::feature::BondKind;
use purr::feature::BracketSymbol;
use purr::feature::Configuration;
use purr::feature::Element;
use purr::feature::VirtualHydrogen;
use purr::walk::Follower;

/// A L-α amino-acid.
#[derive(Clone, Copy, Debug)]
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
}

impl AminoAcid {
    /// Perform a walk on the atoms and bonds of the amino acid.
    ///
    /// The follower must have been initialized with a head. It will finish its
    /// walk on the β carbon, without visiting the atoms part of the peptidic
    /// bond.
    ///
    pub fn visit<F: Follower>(&self, follower: &mut F) {
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

        match self {
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
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH1);
                follower.join(BondKind::Elided, purr::feature::Rnum::try_from(1).unwrap());
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
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.pop(7);
            }

            AminoAcid::Tyr => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.pop(7);
            }

            AminoAcid::Cys => {
                // alpha carbon
                follower.extend(BondKind::Elided, CARBON_TH2);
                // residue
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::S));
                follower.pop(2);
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
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
                follower.pop(1);
                follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::O));
                follower.pop(2);
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
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
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
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::N));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R2);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R1);
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.extend(BondKind::Elided, AtomKind::Aromatic(Aromatic::C));
                follower.join(BondKind::Elided, purr::feature::Rnum::R2);
                follower.pop(10);
            }
        }

        // beta carbon
        follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::C));
    }
}

/// Perform a walk on the atoms and bonds of the protein.
pub fn visit<'aa, S, F>(sequence: S, follower: &mut F)
where
    S: IntoIterator<Item = &'aa AminoAcid>,
    F: Follower,
{
    // visit every amino acid one by one
    let mut aa_iter = sequence.into_iter();
    if let Some(aa) = aa_iter.next() {
        // first amino acid: create a the N of the primary amine and visit residue.
        follower.root(AtomKind::Aliphatic(Aliphatic::N));
        aa.visit(follower);
        // add the carboxy group to the β carbon.
        follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
        follower.pop(1);
        // keep visiting following amino acids.
        while let Some(aa) = aa_iter.next() {
            // next amino acid: create the N atom of the carboxamide and visit residue.
            follower.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
            aa.visit(follower);
            // add the carboxy group to the β carbon.
            follower.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
            follower.pop(1);
        }
        // final amino acid: create the O atom of the carboxylic acid.
        follower.extend(BondKind::Single, AtomKind::Aliphatic(Aliphatic::O));
    }
}

/// Create a SMILES string for the given amino-acid sequence.
pub fn smiles<'aa, S>(sequence: S) -> String
where
    S: IntoIterator<Item = &'aa AminoAcid>,
{
    let mut writer = purr::write::Writer::new();
    visit(sequence, &mut writer);
    writer.write()
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn empty() {
        let s = smiles(&[]);
        assert_eq!(s, "");
    }
}
