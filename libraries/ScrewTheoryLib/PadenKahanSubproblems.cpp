// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#include "ScrewTheoryIkSubproblems.hpp"

#include <cmath>
#include <iostream>

#include "ScrewTheoryTools.hpp"

using namespace roboticslab;

// -----------------------------------------------------------------------------

PadenKahanOne::PadenKahanOne(const MatrixExponential & _exp, const KDL::Vector & _p)
    : exp(_exp),
      p(_p),
      axisPow(vectorPow2(exp.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PadenKahanOne::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - exp.getOrigin();
    KDL::Vector v = k - exp.getOrigin();

    KDL::Vector u_w = axisPow * u;
    KDL::Vector v_w = axisPow * v;

    KDL::Vector u_p = u - u_w;
    KDL::Vector v_p = v - v_w;

    double theta = reference[0];

    if (!KDL::Equal(u_p.Norm(), 0.0) && !KDL::Equal(v_p.Norm(), 0.0))
    {
        theta = std::atan2(KDL::dot(exp.getAxis(), u_p * v_p), KDL::dot(u_p, v_p));
    }

    solutions = {{normalizeAngle(theta)}};

    return KDL::Equal(u_w, v_w) && KDL::Equal(u_p.Norm(), v_p.Norm());
}

// -----------------------------------------------------------------------------

PadenKahanTwo::PadenKahanTwo(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const KDL::Vector & _p, const KDL::Vector & _r)
  : exp1(_exp1),
    exp2(_exp2),
    p(_p),
    r(_r),
    axesCross(exp1.getAxis() * exp2.getAxis()),
    axisPow1(vectorPow2(exp1.getAxis())),
    axisPow2(vectorPow2(exp2.getAxis())),
    axesDot(KDL::dot(exp1.getAxis(), exp2.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PadenKahanTwo::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - r;
    KDL::Vector v = k - r;

    KDL::Vector u_p = u - axisPow2 * u;
    KDL::Vector v_p = v - axisPow1 * v;

    double axis1dot = KDL::dot(exp1.getAxis(), v);
    double axis2dot = KDL::dot(exp2.getAxis(), u);
    double den = std::pow(axesDot, 2) - 1;

    double alpha = (axesDot * axis2dot - axis1dot) / den;
    double beta = (axesDot * axis1dot - axis2dot) / den;

    KDL::Vector term1 = r + alpha * exp1.getAxis() + beta * exp2.getAxis();

    double gamma2 = (std::pow(u.Norm(), 2) - std::pow(alpha, 2) - std::pow(beta, 2) - 2 * alpha * beta * axesDot) / std::pow(axesCross.Norm(), 2);

    bool gamma2_zero = KDL::Equal(gamma2, 0.0);

    bool ret;

    if (!gamma2_zero && gamma2 > 0.0)
    {
        std::cout << "entra en el primer if?\n";
        double gamma = std::sqrt(gamma2);
        KDL::Vector term2 = gamma * axesCross;

        KDL::Vector d = term1 + term2;
        KDL::Vector c = term1 - term2;

        KDL::Vector m = c - r;
        KDL::Vector n = d - r;

        KDL::Vector m1_p = m - axisPow1 * m;
        KDL::Vector m2_p = m - axisPow2 * m;

        KDL::Vector n1_p = n - axisPow1 * n;
        KDL::Vector n2_p = n - axisPow2 * n;
/*
        std::cout << "exp1.getAxis = (" << exp1.getAxis().x() << ", " << exp1.getAxis().y() << ", " << exp1.getAxis().z() << ")\n";
        std::cout << "exp2.getAxis = (" << exp2.getAxis().x() << ", " << exp2.getAxis().y() << ", " << exp2.getAxis().z() << ")\n";
        std::cout << "axesCross = (" << axesCross.x() << ", " << axesCross.y() << ", " << axesCross.z() << ")\n";    
        //std::cout << "axesCross inverted = (" << axesCross_inverted.x() << ", " << axesCross_inverted.y() << ", " << axesCross_inverted.z() << ")\n";    
        //std::cout << "newAxis = (" << newAxes.x() << ", " << newAxes.y() << ", " << newAxes.z() << ")\n";    
        std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
        std::cout << "f - exp2.getOrigin() = (" << (f - exp2.getOrigin()).x() << ", " << (f - exp2.getOrigin()).y() << ", " << (f - exp2.getOrigin()).z() << ")\n";
        std::cout << "exp1.origin = (" << exp1.getOrigin().x() << ", " << exp1.getOrigin().y() << ", " << exp1.getOrigin().z() << ")\n";
        std::cout << "exp2.origin = (" << exp2.getOrigin().x() << ", " << exp2.getOrigin().y() << ", " << exp2.getOrigin().z() << ")\n";    
        //std::cout << "exp4.axis = (" << exp4.getAxis().x() << ", " << exp4.getAxis().y() << ", " << exp4.getAxis().z() << ")\n";    
        //std::cout << "o1 = (" << o1.x() << ", " << o1.y() << ", " << o1.z() << ")\n";
        //std::cout << "o2 = (" << o2.x() << ", " << o2.y() << ", " << o2.z() << ")\n";
       // std::cout << "r4 = (" << r4.x() << ", " << r4.y() << ", " << r4.z() << ")\n";
        std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
        //std::cout << "c2 = (" << c2.x() << ", " << c2.y() << ", " << c2.z() << ")\n";
        //std::cout << "d2 = (" << d2.x() << ", " << d2.y() << ", " << d2.z() << ")\n";
        //std::cout << "pg3_sol_1 = " << pg3_1_sols[0][0] << "\n";
        //std::cout << "pg3_sol_2 = " << pg3_1_sols[1][0] << "\n";
        std::cout << "c = (" << c.x() << ", " << c.y() << ", " << c.z() << ")\n";
        std::cout << "d = (" << d.x() << ", " << d.y() << ", " << d.z() << ")\n";
*/
        double theta1_1 = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        double theta2_1 = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));

        double theta1_2 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        double theta2_2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
/*
        std::cout << "normalize theta1_1 = " << normalizeAngle(theta1_1) <<"\n";
        std::cout << "normalize theta2_1 = " << normalizeAngle(theta2_1) <<"\n";
        std::cout << "normalize theta1_2 = " << normalizeAngle(theta1_2) <<"\n";
        std::cout << "normalize theta2_2 = " << normalizeAngle(theta2_2) <<"\n";
*/
        solutions = {
            {normalizeAngle(theta1_1), normalizeAngle(theta2_1)},
            {normalizeAngle(theta1_2), normalizeAngle(theta2_2)}
        };

        ret = KDL::Equal(m1_p.Norm(), v_p.Norm());
    }
    else
    {
        std::cout << "entra en el else?\n";

        KDL::Vector n = term1 - r;
        KDL::Vector n1_p = n - axisPow1 * n;
        KDL::Vector n2_p = n - axisPow2 * n;

        double theta1 = reference[0];
        double theta2 = reference[1];

        if (!KDL::Equal(v_p.Norm(), 0.0))
        {
            std::cout <<"no entra aqui a que no?\n";
            theta1 = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        }

        if (!KDL::Equal(u_p.Norm(), 0.0))
        {
            std::cout <<"aqui tampoco a que no?\n";
            theta2 = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));
        }

        std::cout << "exp1.getAxis = (" << exp1.getAxis().x() << ", " << exp1.getAxis().y() << ", " << exp1.getAxis().z() << ")\n";
        std::cout << "exp2.getAxis = (" << exp2.getAxis().x() << ", " << exp2.getAxis().y() << ", " << exp2.getAxis().z() << ")\n";
        std::cout << "axesCross = (" << axesCross.x() << ", " << axesCross.y() << ", " << axesCross.z() << ")\n";    
        std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
        std::cout << "f - exp2.getOrigin() = (" << (f - exp2.getOrigin()).x() << ", " << (f - exp2.getOrigin()).y() << ", " << (f - exp2.getOrigin()).z() << ")\n";
        std::cout << "exp1.origin = (" << exp1.getOrigin().x() << ", " << exp1.getOrigin().y() << ", " << exp1.getOrigin().z() << ")\n";
        std::cout << "exp2.origin = (" << exp2.getOrigin().x() << ", " << exp2.getOrigin().y() << ", " << exp2.getOrigin().z() << ")\n";    
        std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
        std::cout << "r = (" << r.x() << ", " << r.y() << ", " << r.z() << ")\n";
        std::cout << "n = (" << n.x() << ", " << n.y() << ", " << n.z() << ")\n";
        std::cout << "n1_p = (" << n1_p.x() << ", " << n1_p.y() << ", " << n1_p.z() << ")\n";
        std::cout << "n2_p = (" << n2_p.x() << ", " << n2_p.y() << ", " << n2_p.z() << ")\n";
        std::cout << "v_p = (" << v_p.x() << ", " << v_p.y() << ", " << v_p.z() << ")\n";
        std::cout << "u_p = (" << u_p.x() << ", " << u_p.y() << ", " << u_p.z() << ")\n";

        double normalized1 = normalizeAngle(theta1);
        double normalized2 = normalizeAngle(theta2);
     
        std::cout << "normalize1 = " << normalized1 <<"\n";
        std::cout << "normalized2 = " << normalized2 <<"\n";

        solutions = {
            {normalized1, normalized2},
            {normalized1, normalized2}
        };

        ret = gamma2_zero && KDL::Equal(n1_p.Norm(), v_p.Norm());
    }

    return ret;
}

// -----------------------------------------------------------------------------

PadenKahanThree::PadenKahanThree(const MatrixExponential & _exp, const KDL::Vector & _p, const KDL::Vector & _k)
    : exp(_exp),
      p(_p),
      k(_k),
      axisPow(vectorPow2(exp.getAxis()))
{}

// -----------------------------------------------------------------------------

bool PadenKahanThree::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector rhsAsVector = rhs * f - k;
    double delta = rhsAsVector.Norm();

    KDL::Vector u = f - exp.getOrigin();
    KDL::Vector v = k - exp.getOrigin();

    KDL::Vector u_p = u - axisPow * u;
    KDL::Vector v_p = v - axisPow * v;

    double alpha = std::atan2(KDL::dot(exp.getAxis(), u_p * v_p), KDL::dot(u_p, v_p));
    double delta_p_2 = std::pow(delta, 2) - std::pow(KDL::dot(exp.getAxis(), f - k), 2);

    double u_p_norm = u_p.Norm();
    double v_p_norm = v_p.Norm();

    bool u_p_norm_zero = KDL::Equal(u_p_norm, 0.0);
    bool v_p_norm_zero = KDL::Equal(v_p_norm, 0.0);

    if (!u_p_norm_zero && !v_p_norm_zero)
    {
        double betaCos = (std::pow(u_p_norm, 2) + std::pow(v_p_norm, 2) - delta_p_2) / (2 * u_p_norm * v_p_norm);
        double betaCosAbs = std::abs(betaCos);
        bool beta_zero_or_pi = KDL::Equal(betaCosAbs, 1.0);

        if (!beta_zero_or_pi && betaCosAbs < 1.0)
        {
            double betaCosCapped = std::max(-1.0, std::min(1.0, betaCos));
            double beta = std::acos(betaCosCapped);

            double theta1 = alpha + beta;
            double theta2 = alpha - beta;

            solutions = {{normalizeAngle(theta1)}, {normalizeAngle(theta2)}};
            return true;
        }
        else
        {
            if (KDL::Equal(betaCos, -1.0))
            {
                alpha += KDL::PI;
            }

            double normalized = normalizeAngle(alpha);
            solutions = {{normalized}, {normalized}};
            return beta_zero_or_pi;
        }
    }
    else
    {
        double ddelta = (f - k).Norm();
        double normalized = normalizeAngle(reference[0]);
        solutions = {{normalized}, {normalized}};
        return KDL::Equal(delta, ddelta);
    }
}

// -----------------------------------------------------------------------------
