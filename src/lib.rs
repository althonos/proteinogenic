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
#[derive(Clone, Copy)]
pub enum AminoAcid {
    Arg,
    His,
    Lys,
    Asp,
    Glu,
    Ser,
    Thr,
    Asn,
    Gln,
    Gly,
    Pro,
    Cys,
    Sec,
    Ala,
    Val,
    Ile,
    Leu,
    Met,
    Phe,
    Tyr,
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

/// Create a SMILES string for the given amino-acid sequence.
pub fn smiles<'aa, S: IntoIterator<Item = &'aa AminoAcid>>(sequence: S) -> String {
    // create a `Writer` to generate the SMILES string.
    let mut writer = purr::write::Writer::new();
    // visit every amino acid one by one
    let mut aa_iter = sequence.into_iter();
    if let Some(aa) = aa_iter.next() {
        // first amino acid: create a the N of the primary amine and visit residue.
        writer.root(AtomKind::Aliphatic(Aliphatic::N));
        aa.visit(&mut writer);
        // add the carboxy group to the β carbon.
        writer.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
        writer.pop(1);
        // keep visiting following amino acids.
        while let Some(aa) = aa_iter.next() {
            // next amino acid: create the N atom of the carboxamide and visit residue.
            writer.extend(BondKind::Elided, AtomKind::Aliphatic(Aliphatic::N));
            aa.visit(&mut writer);
            // add the carboxy group to the β carbon.
            writer.extend(BondKind::Double, AtomKind::Aliphatic(Aliphatic::O));
            writer.pop(1);
        }
        // final amino acid: create the O atom of the carboxylic acid.
        writer.extend(BondKind::Single, AtomKind::Aliphatic(Aliphatic::O));
    }
    // generate final string
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
