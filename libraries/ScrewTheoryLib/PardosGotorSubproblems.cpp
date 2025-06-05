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

PardosGotorEight::PardosGotorEight(const MatrixExponential & _exp1, const MatrixExponential & _exp2, const MatrixExponential & _exp3, const KDL::Vector & _p)
    : exp1(_exp1),
      exp2(_exp2),
      exp3(_exp3),
      p(_p),
      n(computeNormal(exp1, exp2)),
      axisPow(vectorPow2(exp1.getAxis())) // same as exp2.getAxis() and exp3
{}

// -----------------------------------------------------------------------------

bool PardosGotorEight::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector u3 = f - exp3.getOrigin();
    KDL::Vector o3p = exp3.getOrigin() + axisPow * u3;  
    KDL::Vector o3k = rhs * o3p; //ESTÁ BIEN ASÍ????????
    

//------------PONIENDO ESTO ME AHORRARÍA VOLVER A HACER DE NUEVO EL PG4 NO???????-----------------------------------
    //PardosGotorFour pg4(exp1, exp2, o3);
    //pg4.solve(rhs, IDENTIDAD?, reference, solutions);


    //pg4

    KDL::Vector u = o3p - exp2.getOrigin();
    KDL::Vector v = o3k - exp1.getOrigin();

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

    double theta1_co;
    double theta2_oc;
    double theta1_do;
    double theta2_od;

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

        theta1_co = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p), KDL::dot(m1_p, v_p));
        theta2_oc = std::atan2(KDL::dot(exp2.getAxis(), u_p * m2_p), KDL::dot(u_p, m2_p));

        theta1_do = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p), KDL::dot(n1_p, v_p));
        theta2_od = std::atan2(KDL::dot(exp2.getAxis(), u_p * n2_p), KDL::dot(u_p, n2_p));

        /*
        solutions = {
            {normalizeAngle(theta1_co), normalizeAngle(theta2_oc)},
            {normalizeAngle(theta1_do), normalizeAngle(theta2_od)}
        };
        */

        //return samePlane && KDL::Equal(m1_p.Norm(), v_p_norm);//SE CAMBIA POR EL IF DE ABAJO. ESTÁ BIEN?

        if(!(samePlane && KDL::Equal(m1_p.Norm(), v_p_norm))) return false; //HACE LO QUE QUIERO???????????
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

        theta1_co = theta1;
        theta2_oc = theta2;

        theta1_do = theta1_co;
        theta2_od = theta2_oc;

        /*
        double normalized1 = normalizeAngle(theta1);
        double normalized2 = normalizeAngle(theta2);

        solutions = {
            {normalized1, normalized2},
            {normalized1, normalized2}
        };
        */

        //return samePlane && c_zero; //SE CAMBIA POR EL IF DE ABAJO. ESTÁ BIEN?    

        if(!(samePlane && c_zero)) return false; //HACE LO QUE QUIERO???????????
    }

    //pk1

    //MatrixExponential exp3k = joint2twist(w3, o3k, 'rot'); //hay  que usar algo así?
    KDL::Vector pk = o3k + (f - o3p);
    KDL::Vector k = rhs/* * o3p*/ * pk; //BIEN ASI?

    KDL::Vector u = pk - exp3.getOrigin();
    KDL::Vector v = k - exp3.getOrigin();

    KDL::Vector u_w = axisPow * u;
    KDL::Vector v_w = axisPow * v;

    KDL::Vector u_p = u - u_w;
    KDL::Vector v_p = v - v_w;

    double theta = reference[0];

    if (!KDL::Equal(u_p.Norm(), 0.0) && !KDL::Equal(v_p.Norm(), 0.0))
    {
        theta = std::atan2(KDL::dot(exp3.getAxis(), u_p * v_p), KDL::dot(u_p, v_p));
    }

    double theta123 = theta;
    double theta3_1 = theta123 - theta1_co - theta2_oc;
    double theta3_2 = theta123 - theta1_do - theta2_od;

    solutions = {
        {normalizeAngle(theta1_co), normalizeAngle(theta2_oc), normalizeAngle(theta3_1)},
        {normalizeAngle(theta1_do), normalizeAngle(theta2_od), normalizeAngle(theta3_2)}
    };

    return KDL::Equal(u_w, v_w) && KDL::Equal(u_p.Norm(), v_p.Norm());

} 