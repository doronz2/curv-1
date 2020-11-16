	use super::pairing_bls12_381::PairingBls;
	use super::pairing_bls12_381::PAIRING;
	use crate::elliptic::curves::bls12_381::g1;
	use crate::elliptic::curves::bls12_381::g2;
	use crate::BigInt;
	use crate::cryptographic_primitives::hashing::hash_sha256;
	use crate::cryptographic_primitives::hashing::traits::Hash;
	use crate::elliptic::curves::traits::ECPoint;
	use crate::elliptic::curves::traits::ECScalar;
	use serde::export::fmt::Debug;
	use crate::cryptographic_primitives::secret_sharing::feldman_vss;
	use crate::cryptographic_primitives::secret_sharing::feldman_vss::VerifiableSS;
    type GE1 = g1::GE;
    type GE2 = g2::GE;


	//step a in creating phase Ni send commitments to Nj
	pub fn commiting_to_two_polynomials<P:ECPoint>(a_i0: P::Scalar, l: usize, t: usize)-> (VerifiableSS<P>,Vec<P>,Vec<P>){
		let (sss, shares) =
			feldman_vss::VerifiableSS::<P>::share(t,l,a_i0 )  ;
		let commitments_a = sss.commitments;
		let (sss, shares_prime) =
			feldman_vss::VerifiableSS::<P>::share(t,l,P::Scalar::new_random())  ;
		let commitments_b = sss.commitments;
		let commitment_c = commitments_a.iter().zip(commitments_b.iter()).
			map(|comm_a_i,comm_b_i| comm_a_i*comm_b_i).collect();
		(
			VerifiableSS{
				parameters: feldman_vss::ShamirSecretSharing{
					threshold: t,
					share_count: l,
				},
				commitments: commitment_c,
			},
			shares,
			shares_prime
		)
	}

	pub struct Party<P>{
		shares: Vec<P>,
		shares_prime: Vec<P>,
		commitments: Vec<P>,
        index: usize
	}

	impl<P: ECPoint> Party<P> {
        //step a in creating phase Ni send commitments to Nj
        pub fn commiting_to_two_polynomials_to_self<P: ECPoint>(&self, a_i0: P::Scalar, l: usize, t: usize) -> Self {
            let (sss, shares) =
                feldman_vss::VerifiableSS::<P>::share(t, l, a_i0);
            let commitments_a = sss.commitments;
            let (sss, shares_prime) =
                feldman_vss::VerifiableSS::<P>::share(t, l, P::Scalar::new_random());
            let commitments_b = sss.commitments;
            let commitment_c = commitments_a.iter().zip(commitments_b.iter()).
                map(|comm_a_i, comm_b_i| comm_a_i * comm_b_i).collect();
            Self {
                commitments: commitment_c,
                shares,
                shares_prime
            }
        }
        //step b in creating phase
        pub fn validating_commitments(commitment: Vec<P>, shares_j: &[P; 2], index: &u32) -> bool {
            let g: ECPoint = P::generator();
            let h: ECPoint = P::base_point2();
            //computing g^s_ij*h^s'_ij
            let commitment_from_eval: P =
                g.scalar_mul(shares_j[&0]).add_point(h.scalar_mul(shares_j[&1]));
            let mut commitment_iter = commitment.iter();
            let head= commitment_iter.next().unwrap();
            let commitment_from_comms = commitment_iter
                .enumerate()
                .fold(head, |acc, (j,comm)|{
                acc.add_point(comm.scalar_mul( P::Scalar::from(&BigInt::from(index).pow(j as u32))))
            });
            commitment_from_eval == commitment_from_comms
        }
    }
