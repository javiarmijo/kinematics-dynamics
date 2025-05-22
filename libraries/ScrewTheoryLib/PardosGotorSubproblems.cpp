// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#include "ScrewTheoryIkSubproblems.hpp"

#include "ScrewTheoryTools.hpp"

#include <array>

#include <iostream>


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

PardosGotorSix::PardosGotorSix(const MatrixExponential & _exp1, const MatrixExponential & _exp2,
                               const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      p(_p),
      axisPow1(vectorPow2(exp1.getAxis())),
      axisPow2(vectorPow2(exp2.getAxis()))
{}

// -----------------------------------------------------------------------------
bool PardosGotorSix::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform,
                           const JointConfig & reference, Solutions & solutions) const
{
    std::cout <<"entra al solve\n";

    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector r1 = exp1.getOrigin();
    KDL::Vector r2 = exp2.getOrigin();
    KDL::Vector omega1 = exp1.getAxis();
    KDL::Vector omega2 = exp2.getAxis();

    KDL::Vector u = f - r2;
    KDL::Vector v = k - r1;

    KDL::Vector u_p = u - axisPow2 * u;
    KDL::Vector v_p = v - axisPow1 * v;

    KDL::Vector o2 = r2 + axisPow2 * u;
    KDL::Vector o1 = r1 + axisPow1 * v;

    KDL::Vector v3 = omega1 * omega2;

    // r3: intersection line between rotation planes
    const double denom = 1.0 - KDL::dot(omega1, omega2);
    KDL::Vector term1 = omega1 * (KDL::dot(omega1, o1) - KDL::dot(omega2, o2) * KDL::dot(omega1, omega2));
    KDL::Vector term2 = omega2 * (KDL::dot(omega2, o2) - KDL::dot(omega1, o1) * KDL::dot(omega1, omega2));
    KDL::Vector r3 = (term1 + term2) * (1.0 / denom);

    // θ3 for both projections (translation magnitudes)

    // for o2
    const auto diff2 = o2 - r3;
    const double diff2_norm2 = diff2.Norm() * diff2.Norm();
    const double dot2 = KDL::dot(v3, diff2);
    double delta2_2 = dot2 * dot2 - diff2_norm2 + u_p.Norm() * u_p.Norm();
    //if (delta2_2 < 0.0) return false;
    const double delta2 = std::sqrt(delta2_2);

    const double t2a = dot2 + delta2;
    const double t2b = dot2 - delta2;
    const auto cand2a = r3 + t2a * v3;
    const auto cand2b = r3 + t2b * v3;

    double θ3_c2, θ3_d2;
    if ((cand2a - o2).Norm() < (cand2b - o2).Norm()) {
        θ3_c2 = t2a;
        θ3_d2 = t2b;
    } else {
        θ3_c2 = t2b;
        θ3_d2 = t2a;
    }

    // for o1
    const auto diff1 = o1 - r3;
    const double diff1_norm2 = diff1.Norm() * diff1.Norm();
    const double dot1 = KDL::dot(v3, diff1);
    double delta1_2 = dot1 * dot1 - diff1_norm2 + v_p.Norm() * v_p.Norm();
    //if (delta1_2 < 0.0) return false;
    const double delta1 = std::sqrt(delta1_2);

    const double t1a = dot1 + delta1;
    const double t1b = dot1 - delta1;
    const auto cand1a = r3 + t1a * v3;
    const auto cand1b = r3 + t1b * v3;

    double θ3_c1, θ3_d1;
    if ((cand1a - o1).Norm() < (cand1b - o1).Norm()) {
        θ3_c1 = t1a;
        θ3_d1 = t1b;
    } else {
        θ3_c1 = t1b;
        θ3_d1 = t1a;
    }

    const auto c2 = r3 + θ3_c2 * v3;
    const auto d2 = r3 + θ3_d2 * v3;
    const auto c1 = r3 + θ3_c1 * v3;
    const auto d1 = r3 + θ3_d1 * v3;

    bool delta1_2_zero = KDL::Equal(delta1_2, 0.0);
    bool delta2_2_zero = KDL::Equal(delta2_2, 0.0);

    bool ret;

    if(!delta1_2_zero && delta1_2 > 0.0 && !delta2_2_zero && delta2_2 > 0.0)
    {
        std::cout <<"entra al if\n";
        std::cout <<"delta1_2:";
        std::cout <<delta1_2;
        std::cout <<"\n";
        std::cout <<"delta2_2:";
        std::cout <<delta2_2;
        std::cout <<"\n";

        const auto m2 = c2 - r2;
        const auto n2 = d2 - r2;
        const auto m1 = c1 - r1;
        const auto n1 = d1 - r1;

        const auto mp2 = m2 - axisPow2 * m2;
        const auto np2 = n2 - axisPow2 * n2;
        const auto mp1 = m1 - axisPow1 * m1;
        const auto np1 = n1 - axisPow1 * n1;

        auto pk1Theta = [](const KDL::Vector & omega, const KDL::Vector & a, const KDL::Vector & b) {
            return std::atan2(KDL::dot(omega, a * b), KDL::dot(a, b));
        };

        const double θ2_c = pk1Theta(omega2, u_p, mp2);
        const double θ2_d = pk1Theta(omega2, u_p, np2);
        const double θ1_c = pk1Theta(omega1, mp1, v_p);
        const double θ1_d = pk1Theta(omega1, np1, v_p);

        solutions = {
            { normalizeAngle(θ1_c), normalizeAngle(θ2_c) },
            { normalizeAngle(θ1_d), normalizeAngle(θ2_d) }
        };


        
        ret = KDL::Equal(mp1.Norm(), v_p.Norm());
    }

    else
    {
        std::cout <<"entra al else\n";
                std::cout <<"delta1_2:";
        std::cout <<delta1_2;
        std::cout <<"\n";
        std::cout <<"delta2_2:";
        std::cout <<delta2_2;
        std::cout <<"\n";

        bool valid = true;
        
        // θ3 for both projections (translation magnitudes)

        // for o2
        if (delta2_2 < 0.0 && delta2_2_zero) delta2_2=0.00;
        else if(delta2_2 < 0.0 && !delta2_2_zero)
        {
            delta2_2=0.00;
            valid=false;
        } 
        
        const double delta2 = std::sqrt(delta2_2);
        const double t2a = dot2 + delta2;
        const double t2b = dot2 - delta2;
        const auto cand2a = r3 + t2a * v3;
        const auto cand2b = r3 + t2b * v3;

        double θ3_c2, θ3_d2;
        if ((cand2a - o2).Norm() < (cand2b - o2).Norm()) {
            θ3_c2 = t2a;
            θ3_d2 = t2b;
        } else {
            θ3_c2 = t2b;
            θ3_d2 = t2a;
        }

        // for o1
        if (delta1_2 < 0.0 && delta1_2_zero) delta1_2=0.0;       
        else if(delta1_2 < 0.0 && !delta1_2_zero)
        {
            valid=false;
            delta1_2=0.0;
        } 

        const double delta1 = std::sqrt(delta1_2);  
        const double t1a = dot1 + delta1;
        const double t1b = dot1 - delta1;
        const auto cand1a = r3 + t1a * v3;
        const auto cand1b = r3 + t1b * v3;

        double θ3_c1, θ3_d1;
        if ((cand1a - o1).Norm() < (cand1b - o1).Norm()) {
            θ3_c1 = t1a;
            θ3_d1 = t1b;
        } else {
            θ3_c1 = t1b;
            θ3_d1 = t1a;
        }

        const auto c2 = r3 + θ3_c2 * v3;
        const auto d2 = r3 + θ3_d2 * v3;
        const auto c1 = r3 + θ3_c1 * v3;
        const auto d1 = r3 + θ3_d1 * v3;

        double θ2_c = reference[1];
        double θ2_d = reference[1];
        double θ1_c = reference[0];
        double θ1_d = reference[0];

        //Aproximate solution
        const auto m2 = c2 - r2;
        const auto n2 = d2 - r2;
        const auto m1 = c2 - r1;//cambiado para el check solutions 2
        const auto n1 = d2 - r1;//cambiado para el chekc solutions 2

        const auto mp2 = m2 - axisPow2 * m2;
        const auto np2 = n2 - axisPow2 * n2;
        const auto mp1 = m1 - axisPow1 * m1;
        const auto np1 = n1 - axisPow1 * n1;


        auto pk1Theta = [](const KDL::Vector & omega, const KDL::Vector & a, const KDL::Vector & b) {
            return std::atan2(KDL::dot(omega, a * b), KDL::dot(a, b));
        };

        if (!KDL::Equal(u_p.Norm(), 0.0))
        {
            θ2_c = pk1Theta(omega2, u_p, mp2);
            θ2_d = pk1Theta(omega2, u_p, np2);
        }
       
        if (!KDL::Equal(v_p.Norm(), 0.0))
        {
            θ1_c = pk1Theta(omega1, mp1, v_p);
            θ1_d = pk1Theta(omega1, np1, v_p);            
        }


        solutions = {
            { normalizeAngle(θ1_c), normalizeAngle(θ2_c) },
            { normalizeAngle(θ1_d), normalizeAngle(θ2_d) }
        };



        ret = (delta1_2_zero || delta2_2_zero) && (KDL::Equal(np1.Norm(), v_p.Norm()) && valid);
       

    }
    
    return ret;
}
