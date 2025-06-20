// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#include "ScrewTheoryIkSubproblems.hpp"

#include "ScrewTheoryTools.hpp"

using namespace roboticslab;

// -----------------------------------------------------------------------------

namespace
{
    KDL::Vector computeNormal(const MatrixExponential & exp1, const MatrixExponential & exp2)
    {
        KDL::Vector diff = exp2.getOrigin() - exp1.getOrigin();
        KDL::Vector normal = (exp1.getAxis() * diff) * exp1.getAxis();
        normal.Normalize();
        return vectorPow2(normal) * diff;
    }
}

// -----------------------------------------------------------------------------

PardosGotorOne::PardosGotorOne(const MatrixExponential & _exp, const KDL::Vector & _p)
    : exp(_exp),
      p(_p)
{}

// -----------------------------------------------------------------------------

bool PardosGotorOne::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector diff = k - f;
    double theta = KDL::dot(exp.getAxis(), diff);

    solutions = {{theta}};
    return true;
}

// -----------------------------------------------------------------------------

PardosGotorTwo::PardosGotorTwo(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      crossPr2(exp2.getAxis() * exp1.getAxis()),
      crossPr2Norm(crossPr2.Norm())
{}

// -----------------------------------------------------------------------------

bool PardosGotorTwo::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector crossPr1 = exp2.getAxis() * (f - k);
    double crossPr1Norm = crossPr1.Norm();

    KDL::Vector c;

    if (KDL::dot(crossPr1, crossPr2) >= crossPr1Norm * crossPr2Norm)
    {
        c = k + (crossPr1Norm / crossPr2Norm) * exp1.getAxis();
    }
    else
    {
        c = k - (crossPr1Norm / crossPr2Norm) * exp1.getAxis();
    }

    double theta1 = KDL::dot(exp1.getAxis(), k - c);
    double theta2 = KDL::dot(exp2.getAxis(), c - f);

    solutions = {{theta1, theta2}};

    return true;
}

// -----------------------------------------------------------------------------

PardosGotorThree::PardosGotorThree(const MatrixExponential & _exp, const KDL::Vector & _p, const KDL::Vector & _k)
    : exp(_exp),
      p(_p),
      k(_k)
{}

// -----------------------------------------------------------------------------

bool PardosGotorThree::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector rhsAsVector = rhs * p - k;
    double delta = rhsAsVector.Norm();

    KDL::Vector diff = k - f;

    double dotPr = KDL::dot(exp.getAxis(), diff);
    double sq2 = std::pow(dotPr, 2) - std::pow(diff.Norm(), 2) + std::pow(delta, 2);
    bool sq2_zero = KDL::Equal(sq2, 0.0);

    bool ret;

    if (!sq2_zero && sq2 > 0)
    {
        double sq = std::sqrt(std::abs(sq2));
        solutions = {{dotPr + sq}, {dotPr - sq}};
        ret = true;
    }
    else
    {
        KDL::Vector proy = vectorPow2(exp.getAxis()) * diff;
        double norm = proy.Norm();
        solutions = {{norm}, {norm}};
        ret = sq2_zero;
    }

    return ret;
}

// -----------------------------------------------------------------------------

PardosGotorFour::PardosGotorFour(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      n(computeNormal(exp1, exp2)),
      axisPow(vectorPow2(exp1.getAxis())) // same as exp2.getAxis()
{}

// -----------------------------------------------------------------------------

bool PardosGotorFour::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - exp2.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p = u - axisPow * u;
    KDL::Vector v_p = v - axisPow * v;

    KDL::Vector c1 = exp1.getOrigin() + v - v_p;
    KDL::Vector c2 = exp2.getOrigin() + u - u_p;

    KDL::Vector c_diff = c2 - c1;
    bool samePlane = KDL::Equal(c_diff, n);

    if (!samePlane)
    {
        c_diff = n; // proyection of c_diff onto the perpendicular plane
        c1 = c2 - c_diff; // c1 on the intersecion of axis 1 and the normal plane to both axes
    }

    double c_norm = c_diff.Norm();
    double u_p_norm = u_p.Norm();
    double v_p_norm = v_p.Norm();

    double c_test = u_p_norm + v_p_norm - c_norm;
    bool c_zero = KDL::Equal(c_test, 0.0);

    if (!c_zero && c_test > 0.0 && u_p_norm > 0.0 && v_p_norm > 0.0)
    {
        KDL::Vector omega_a = c_diff / c_norm;
        KDL::Vector omega_h = exp1.getAxis() * omega_a;

        double a = (std::pow(c_norm, 2) - std::pow(u_p_norm, 2) + std::pow(v_p_norm, 2)) / (2 * c_norm);
        double h = std::sqrt(std::abs(std::pow(v_p.Norm(), 2) - std::pow(a, 2)));

        KDL::Vector term1 = c1 + a * omega_a;
        KDL::Vector term2 = h * omega_h;

        KDL::Vector c = term1 + term2;
        KDL::Vector d = term1 - term2;

        KDL::Vector m1 = c - exp1.getOrigin();
        KDL::Vector m2 = c - exp2.getOrigin();

        KDL::Vector n1 = d - exp1.getOrigin();
        KDL::Vector n2 = d - exp2.getOrigin();

        KDL::Vector m1_p = m1 - axisPow * m1;
        KDL::Vector m2_p = m2 - axisPow * m2;

        KDL::Vector n1_p = n1 - axisPow * n1;
        KDL::Vector n2_p = n2 - axisPow * n2;

        double theta1_1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        double theta2_1 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));

        double theta1_2 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        double theta2_2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));

        solutions = {
            {normalizeAngle(theta1_1), normalizeAngle(theta2_1)},
            {normalizeAngle(theta1_2), normalizeAngle(theta2_2)}
        };

        return samePlane && KDL::Equal(m1_p.Norm(), v_p_norm);
    }
    else
    {
        double theta1 = reference[0];
        double theta2 = reference[1];

        if (!KDL::Equal(v_p_norm, 0.0))
        {
            theta1 = std::atan2(KDL::dot(exp1.getAxis(), c_diff * v_p), KDL::dot(c_diff, v_p));
        }

        if (!KDL::Equal(u_p_norm, 0.0))
        {
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * c_diff), KDL::dot(-c_diff, u_p));
        }

        double normalized1 = normalizeAngle(theta1);
        double normalized2 = normalizeAngle(theta2);

        solutions = {
            {normalized1, normalized2},
            {normalized1, normalized2}
        };

        return samePlane && c_zero;
    }
}

// -----------------------------------------------------------------------------

PardosGotorSix::PardosGotorSix(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      axesCross(exp1.getAxis() * exp2.getAxis()),
      axisPow1(vectorPow2(exp1.getAxis())),
      axisPow2(vectorPow2(exp2.getAxis())),
      axesDot(KDL::dot(exp1.getAxis(), exp2.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorSix::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - exp2.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p = u - axisPow2 * u;
    KDL::Vector v_p = v - axisPow1 * v;

    KDL::Vector o2 = exp2.getOrigin() + axisPow2 * u;
    KDL::Vector o1 = exp1.getOrigin() + axisPow1 * v;

    double o2_dot = KDL::dot(exp2.getAxis(), o2);
    double o1_dot = KDL::dot(exp1.getAxis(), o1);

    KDL::Vector r3 = (exp1.getAxis() * (o1_dot - o2_dot * axesDot) + exp2.getAxis() * (o2_dot - o1_dot * axesDot)) / (1 - axesDot);

    MatrixExponential exp3(MatrixExponential::TRANSLATION, axesCross);
    PardosGotorThree pg3_2(exp3, r3, o2);
    PardosGotorThree pg3_1(exp3, r3, o1);

    Solutions pg3_1_sols, pg3_2_sols;

    bool pg3_1_ret = pg3_1.solve(KDL::Frame(v_p - (r3 - o1)), KDL::Frame::Identity(), pg3_1_sols);
    bool pg3_2_ret = pg3_2.solve(KDL::Frame(u_p - (r3 - o2)), KDL::Frame::Identity(), pg3_2_sols);

    KDL::Vector c2 = r3 + pg3_2_sols[0][0] * exp3.getAxis();
    KDL::Vector d2 = r3 + pg3_2_sols[1][0] * exp3.getAxis();

    KDL::Vector c1 = r3 + pg3_1_sols[0][0] * exp3.getAxis();
    KDL::Vector d1 = r3 + pg3_1_sols[1][0] * exp3.getAxis();

    bool ret = pg3_1_ret && pg3_2_ret;

    double theta1, theta2;

    if (c1 == c2)
    {
        KDL::Vector m2 = c2 - exp2.getOrigin();
        KDL::Vector m1 = c1 - exp1.getOrigin();

        KDL::Vector m2_p = m2 - axisPow2 * m2;
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        theta1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));
    }
    else if (d1 == d2)
    {
        KDL::Vector n2 = d2 - exp2.getOrigin();
        KDL::Vector n1 = d1 - exp1.getOrigin();

        KDL::Vector n2_p = n2 - axisPow2 * n2;
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        theta1 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
    }
    else
    {
        // these might be equal if `pg3_2_ret` is false
        KDL::Vector m1 = c2 - exp1.getOrigin();
        KDL::Vector n1 = d2 - exp1.getOrigin();

        if (m1.Norm() <= n1.Norm())
        {
            KDL::Vector m2 = c2 - exp2.getOrigin();

            KDL::Vector m2_p = m2 - axisPow2 * m2;
            KDL::Vector m1_p = m1 - axisPow1 * m1;

            theta1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));
        }
        else
        {
            KDL::Vector n2 = d2 - exp2.getOrigin();

            KDL::Vector n2_p = n2 - axisPow2 * n2;
            KDL::Vector n1_p = n1 - axisPow1 * n1;

            theta1 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
        }

        if(f != k) ret = false;
    }

    solutions = {{normalizeAngle(theta1), normalizeAngle(theta2)}};

    return ret;
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

PardosGotorSeven::PardosGotorSeven(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const MatrixExponential & _exp3, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      exp3(_exp3),
      p(_p),
      n(computeNormal(exp1, exp2)),
      axisPow1(vectorPow2(exp1.getAxis())),
      axisPow2(vectorPow2(exp2.getAxis())),
      axesCross(exp1.getAxis() * exp2.getAxis()),
      axesCross_inverted(exp2.getAxis() * exp1.getAxis()),
      axesDot(KDL::dot(exp1.getAxis(), exp2.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PardosGotorSeven::solve(const KDL::Frame &rhs, const KDL::Frame &pointTransform, const JointConfig &reference, Solutions &solutions) const {
    
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;
    KDL::Vector u = f - exp3.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p2 = u - axisPow2 * u;
    KDL::Vector v_p1 = v - axisPow1 * v;

    KDL::Vector o2 = exp2.getOrigin() + axisPow2 * (f - exp2.getOrigin());
    KDL::Vector o1 = exp1.getOrigin() + axisPow1 * v;

    double o2_dot = KDL::dot(exp2.getAxis(), o2);
    double o1_dot = KDL::dot(exp1.getAxis(), o1);

    KDL::Vector r4 = (exp1.getAxis() * (o1_dot - o2_dot * axesDot) + exp2.getAxis() * (o2_dot - o1_dot * axesDot)) / (1 - axesDot);
    KDL::Vector newAxes = KDL::dot(axesCross, o1 - o2) < 0.0 ? axesCross_inverted : axesCross;

    MatrixExponential exp4(MatrixExponential::TRANSLATION, newAxes);
    PardosGotorThree pg3_1(exp4, r4, o1);

    Solutions pg3_1_sols;
    if (!pg3_1.solve(KDL::Frame(v_p1 - (r4 - o1)), KDL::Frame::Identity(), pg3_1_sols)) return false;

    KDL::Vector c1 = r4 + (pg3_1_sols[0][0] * exp4.getAxis());
    KDL::Vector d1 = r4 + (pg3_1_sols[1][0] * exp4.getAxis());

    double theta_ck = reference[0], theta_dk = reference[1];
    PardosGotorFour pg4(exp2, exp3, f);

    Solutions pg4_c_sols, pg4_d_sols;
    bool pg4_ret_c = pg4.solve(KDL::Frame(c1 - f), KDL::Frame::Identity(), reference, pg4_c_sols);
    bool pg4_ret_d = pg4.solve(KDL::Frame(d1 - f), KDL::Frame::Identity(), reference, pg4_d_sols);

    if ((!pg4_ret_c && !pg4_ret_d) && (KDL::dot(axesCross, o1 - o2) == 0.0))
    {
        c1 = r4 + (-pg3_1_sols[0][0] * exp4.getAxis());
        d1 = r4 + (-pg3_1_sols[1][0] * exp4.getAxis());

        pg4_ret_c = pg4.solve(KDL::Frame(c1 - f), KDL::Frame::Identity(), reference, pg4_c_sols);
        pg4_ret_d = pg4.solve(KDL::Frame(d1 - f), KDL::Frame::Identity(), reference, pg4_d_sols);
    }

    if (pg4_ret_c && pg4_ret_d)
    {
        KDL::Vector m1 = c1 - exp1.getOrigin();
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p1), KDL::dot(m1_p, v_p1));            
        }
    }
    else if (pg4_ret_c)
    {
        KDL::Vector m1 = c1 - exp1.getOrigin();
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_ck = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p1), KDL::dot(m1_p, v_p1));
            theta_dk = theta_ck;            
        }

        pg4_d_sols = pg4_c_sols;
    }
    else if (pg4_ret_d)
    {
        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = theta_dk;            
        }

        pg4_c_sols = pg4_d_sols;
    }
    else
    {
        return false;
    } 

    solutions = {
        {theta_ck, pg4_c_sols[0][0], pg4_c_sols[0][1]},
        {theta_ck, pg4_c_sols[1][0], pg4_c_sols[1][1]},
        {theta_dk, pg4_d_sols[0][0], pg4_d_sols[0][1]},
        {theta_dk, pg4_d_sols[1][0], pg4_d_sols[1][1]}
    };

    return true;
}


// -----------------------------------------------------------------------------