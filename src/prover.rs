use ark_ff::FftField;
use ark_poly::{EvaluationDomain, univariate::DensePolynomial};
use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment, LabeledCommitment, QuerySet};
use ark_std::rand::RngCore;
use ark_ff::to_bytes;
use ark_poly::{UVPolynomial};

use crate::{
    data_structures::{Proof, Statement},
    error::Error,
    rng::FiatShamirRng,
    PROTOCOL_NAME,
};

pub fn prove<
    F: FftField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
    FS: FiatShamirRng,
    R: RngCore,
>(
    ck: &PC::CommitterKey,
    statement: &Statement<F, PC>,
    f_poly_l: &LabeledPolynomial<F, DensePolynomial<F>>,
    f_rand: &PC::Randomness,
    rng: &mut R,
) -> Result<Proof<F, PC>, Error<PC::Error>> {
    // In the rest of protocol that is not described here, the masking polynomial is opened twice. Therefore, the masking polynomial cannot be a constant polynomial.
    let domain = statement.domain;
    let _sum = statement.sum;
    let f: PC::Commitment = statement.f.clone();
    
    let mut fs_rng = FS::initialize(&to_bytes![&PROTOCOL_NAME, statement].unwrap());
    let h_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    let g_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    let s_poly = f_poly_l.polynomial() * (-F::from(1u64));
    
    let s_poly_l = LabeledPolynomial::new("s".into(), s_poly, None, None);
    let h_poly_l = LabeledPolynomial::new("h".into(), h_poly, None, None);
    let g_poly_l = LabeledPolynomial::new("g".into(), g_poly, Some(domain.size() - 2), None);
    let (s_commitment, s_rand) = PC::commit(ck, &[s_poly_l.clone()], Some(rng)).unwrap();
    let (h_commitment, h_rand) = PC::commit(ck, &[h_poly_l.clone()], Some(rng)).unwrap();
    let (g_commitment, g_rand) = PC::commit(ck, &[g_poly_l.clone()], Some(rng)).unwrap();
    

    let s: PC::Commitment = s_commitment[0].commitment().clone();
    let h: PC::Commitment = h_commitment[0].commitment().clone();
    let g: PC::Commitment = g_commitment[0].commitment().clone();

    fs_rng.absorb(&to_bytes![s, h, g].unwrap());

    let xi = F::rand(&mut fs_rng);
    let opening_challenge = F::rand(&mut fs_rng);

    let f_opening = f_poly_l.evaluate(&xi);
    let s_opening = s_poly_l.evaluate(&xi);
    let g_opening = g_poly_l.evaluate(&xi);
    let h_opening = h_poly_l.evaluate(&xi);

    let f_l = LabeledCommitment::new("s".into(), f.clone(), None);
    let s_l = LabeledCommitment::new("s".into(), s.clone(), None);
    let h_l = LabeledCommitment::new("h".into(), h.clone(), None);
    let g_l = LabeledCommitment::new("g".into(), g.clone(), Some(domain.size()-2));
    
    let point_label = String::from("xi");
    let query_set = QuerySet::from([
        ("f".into(), (point_label.clone(), xi)),
        ("h".into(), (point_label.clone(), xi)),
        ("g".into(), (point_label.clone(), xi)),
        ("s".into(), (point_label, xi)),
    ]);

    let pc_proof = PC::batch_open(
        ck,
        [f_poly_l, &s_poly_l, &g_poly_l, &h_poly_l],
        &[f_l, s_l, h_l, g_l],
        &query_set,
        opening_challenge,
        [f_rand, &s_rand[0], &g_rand[0], &h_rand[0]],
        Some(rng)
    ).unwrap();
    return Ok(Proof {
        f_opening,
        s,
        s_opening,
        h,
        h_opening,
        g,
        g_opening,
        pc_proof,
    });
}

// s(xi) + f(xi) = xi * g(xi) + h(xi) * pi(xi-d_i) + sum * d_s^-1
// s(xi) + f(xi) = xi * g(xi) + h(xi) * pi(xi-d_i) + sum * (predictable random number)
// s(xi) + f(xi) = xi * g(xi) + h(xi) * pi(xi-d_i) + y
// s(x) + f(x) = x * g(x) + h(x) * pi(x-d_i) + y

// f(x) = beta + x * g(x) + h(x) * pi(x-d_i)
// s(x) = y - beta