// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#include "ScrewTheoryIkSubproblems.hpp"

#include "ScrewTheoryTools.hpp"

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
    std::cout << "sq2 = " << sq2 << "\n";
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
        std::cout <<"pg4 1\n";
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

        std::cout << "equal = " <<KDL::Equal(m1_p.Norm(), v_p_norm)<<"\n";
        return samePlane && KDL::Equal(m1_p.Norm(), v_p_norm);
    }
    else
    {
        std::cout <<"pg4 2\n";
        double theta1 = reference[0];
        double theta2 = reference[1];

        if (!KDL::Equal(v_p_norm, 0.0))
        {
            std::cout <<"si no?\n";
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

bool PardosGotorSeven::solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const
{
    KDL::Vector f = pointTransform * p;
    KDL::Vector k = rhs * p;

    KDL::Vector u = f - exp3.getOrigin();
    KDL::Vector v = k - exp1.getOrigin();

    KDL::Vector u_p2 = (f - exp2.getOrigin()) - axisPow2 * (f - exp2.getOrigin());
    KDL::Vector v_p1 = v - axisPow1 * v;

    KDL::Vector o2 = exp2.getOrigin() + axisPow2 * (f - exp2.getOrigin());
    //KDL::Vector o3 = exp3.getOrigin() + axisPow3 * u;
    KDL::Vector o1 = exp1.getOrigin() + axisPow1 * v;

    double o2_dot = KDL::dot(exp2.getAxis(), o2);
    //double o3_dot = KDL::dot(exp3.getAxis(), o3);
    double o1_dot = KDL::dot(exp1.getAxis(), o1);

    KDL::Vector r4 = (exp1.getAxis() * (o1_dot - o2_dot * axesDot) + exp2.getAxis() * (o2_dot - o1_dot * axesDot)) / (1 - axesDot);
    //KDL::Vector r4 = (exp1.getAxis() * (o1_dot - o3_dot * axesDot) + exp3.getAxis() * (o3_dot - o1_dot * axesDot)) / (1 - axesDot);

    KDL::Vector newAxes = axesCross;
    KDL::Vector dir = o1 - o2;

    std::cout << "dot(axesCross, dir) = " << dot(axesCross, dir) << "\n";

    if (KDL::dot(axesCross, dir) < 0.0)
    {
        std::cout <<"hola\n";
        newAxes = axesCross_inverted;  // invertir si está en sentido opuesto al giro real
    }
        
/**/MatrixExponential exp4(MatrixExponential::TRANSLATION, newAxes/*MAAAAL --- SI PONGO NEW AXES DA SEGFAULT*/);//A VECES COGE SENTIDO CONTRARIO. CORREGIR
    PardosGotorThree pg3_1(exp4, r4, o1);
    //PardosGotorThree pg3_2(exp4, r4, o2);

    Solutions pg3_1_sols, pg3_2_sols;

    bool pg3_1_ret = pg3_1.solve(KDL::Frame(v_p1 - (r4 - o1)), KDL::Frame::Identity(), pg3_1_sols);
    //bool pg3_2_ret = pg3_2.solve(KDL::Frame(u_p2 - (r4 - o2)), KDL::Frame::Identity(), pg3_2_sols);

    bool ret = pg3_1_ret; 

    if(!ret) return false;

  //KDL::Vector c2 = r4 - pg3_1_sols[0][0] * exp4.getAxis();
  //KDL::Vector d2 = r4 - pg3_1_sols[1][0] * exp4.getAxis();

/**/KDL::Vector c1 = r4 + pg3_1_sols[0][0] * exp4.getAxis();//EL AXIS A VECES SALE EN SENTIDO CONTRARIO AL ESPERADO Y POR ESO EL ERROR
    KDL::Vector d1 = r4 + pg3_1_sols[1][0] * exp4.getAxis();

    std::cout << "exp1.getAxis = (" << exp1.getAxis().x() << ", " << exp1.getAxis().y() << ", " << exp1.getAxis().z() << ")\n";
    std::cout << "exp2.getAxis = (" << exp2.getAxis().x() << ", " << exp2.getAxis().y() << ", " << exp2.getAxis().z() << ")\n";
    std::cout << "axesCross = (" << axesCross.x() << ", " << axesCross.y() << ", " << axesCross.z() << ")\n";    
    std::cout << "axesCross inverted = (" << axesCross_inverted.x() << ", " << axesCross_inverted.y() << ", " << axesCross_inverted.z() << ")\n";    
    std::cout << "newAxis = (" << newAxes.x() << ", " << newAxes.y() << ", " << newAxes.z() << ")\n";    
    std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
    std::cout << "f - exp2.getOrigin() = (" << (f - exp2.getOrigin()).x() << ", " << (f - exp2.getOrigin()).y() << ", " << (f - exp2.getOrigin()).z() << ")\n";
    std::cout << "exp1.origin = (" << exp1.getOrigin().x() << ", " << exp1.getOrigin().y() << ", " << exp1.getOrigin().z() << ")\n";
    std::cout << "exp2.origin = (" << exp2.getOrigin().x() << ", " << exp2.getOrigin().y() << ", " << exp2.getOrigin().z() << ")\n";    
    std::cout << "exp4.axis = (" << exp4.getAxis().x() << ", " << exp4.getAxis().y() << ", " << exp4.getAxis().z() << ")\n";    
    std::cout << "o1 = (" << o1.x() << ", " << o1.y() << ", " << o1.z() << ")\n";
    std::cout << "o2 = (" << o2.x() << ", " << o2.y() << ", " << o2.z() << ")\n";
    std::cout << "r4 = (" << r4.x() << ", " << r4.y() << ", " << r4.z() << ")\n";
    std::cout << "k = (" << k.x() << ", " << k.y() << ", " << k.z() << ")\n";
    //std::cout << "c2 = (" << c2.x() << ", " << c2.y() << ", " << c2.z() << ")\n";
    //std::cout << "d2 = (" << d2.x() << ", " << d2.y() << ", " << d2.z() << ")\n";
    std::cout << "pg3_sol_1 = " << pg3_1_sols[0][0] << "\n";
    std::cout << "pg3_sol_2 = " << pg3_1_sols[1][0] << "\n";
    std::cout << "c1 = (" << c1.x() << ", " << c1.y() << ", " << c1.z() << ")\n";
    std::cout << "d1 = (" << d1.x() << ", " << d1.y() << ", " << d1.z() << ")\n";

    //double theta_ck, theta_dk;

    double theta_dk = reference[0];
    double theta_ck = reference[1];

    PardosGotorFour pg4(exp2, exp3, f);

    Solutions pg4_c_sols, pg4_d_sols;
    bool pg4_ret_c, pg4_ret_d;

    pg4_ret_c = pg4.solve(KDL::Frame(c1 - f), KDL::Frame::Identity(), reference, pg4_c_sols);

    pg4_ret_d = pg4.solve(KDL::Frame(d1 - f), KDL::Frame::Identity(), reference, pg4_d_sols);

    std::cout << "pg4_ret_c = " << pg4_ret_c <<" | pg4_ret_d = " << pg4_ret_d << "\n";

    if (pg4_ret_c && pg4_ret_d)
    {
        /*
        PadenKahanOne pk1c(exp1, c1);
        PadenKahanOne pk1d(exp1, d1);

        Solutions pk1_sol_c, pk1_sol_d;
        bool pk1_ret_c, pk1_ret_d;

        pk1_ret_c = pk1c.solve(KDL::Frame(k - c1), KDL::Frame::Identity(), reference, pk1_sol_c);
        pk1_ret_d = pk1d.solve(KDL::Frame(k - d1), KDL::Frame::Identity(), reference, pk1_sol_d);

        std::cout << " pk1_ret_d = " << pk1_ret_d << "\n";

        if(!pk1_ret_d && d1!=k) return false;
        else if(!pk1_ret_c && c1!=k) return false;

        theta_ck = pk1_sol_c[0][0];
        theta_dk = pk1_sol_d[0][0];
        */
        
        std::cout << "no entra no?\n";
        KDL::Vector m1 = c1 - exp1.getOrigin();
        KDL::Vector m1_p = m1 - axisPow1 * m1;

        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        std::cout << "v = (" << v.x() << ", " << v.y() << ", " << v.z() << ")\n";
        std::cout << "v_p1 = (" << v_p1.x() << ", " << v_p1.y() << ", " << v_p1.z() << ")\n";
        std::cout << "n1 = (" << n1.x() << ", " << n1.y() << ", " << n1.z() << ")\n";
        std::cout << "n1_p = (" << n1_p.x() << ", " << n1_p.y() << ", " << n1_p.z() << ")\n";
        std::cout << "m1 = (" << m1.x() << ", " << m1.y() << ", " << m1.z() << ")\n";
        std::cout << "m1_p = (" << m1_p.x() << ", " << m1_p.y() << ", " << m1_p.z() << ")\n";

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = std::atan2(KDL::dot(exp1.getAxis(), m1_p * v_p1), KDL::dot(m1_p, v_p1));            
        }

        std::cout << "theta_ck = " << theta_ck <<"\n";
        std::cout << "theta_dk = " << theta_dk <<"\n";
        
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
        /*
        PadenKahanOne pk1(exp1, d1);

        Solutions pk1_sol_d;
        bool pk1_ret_d;

        pk1_ret_d = pk1.solve(KDL::Frame(k - d1), KDL::Frame::Identity(), reference, pk1_sol_d);

        std::cout << " pk1_ret_d = " << pk1_ret_d << "\n";

        if(!pk1_ret_d) return false;
        */
        ///*
        KDL::Vector n1 = d1 - exp1.getOrigin();
        KDL::Vector n1_p = n1 - axisPow1 * n1;

        if (!KDL::Equal(v_p1.Norm(), 0.0))
        {
            theta_dk = std::atan2(KDL::dot(exp1.getAxis(), n1_p * v_p1), KDL::dot(n1_p, v_p1));
            theta_ck = theta_dk;            
        }
        //*/

        /*
        theta_dk = pk1_sol_d[0][0];
        theta_ck = pk1_sol_d[0][0];
        */
        pg4_c_sols = pg4_d_sols;
    }
    else
    {
        std::cout << "problema aquí?\n";
        return false;
    } 

        std::cout << "theta_ck = " << theta_ck <<"\n";
        std::cout << "theta_dk = " << theta_dk <<"\n";

    solutions = {
        {theta_ck, pg4_c_sols[0][0], pg4_c_sols[0][1]},    // las soluciones 1 y 3 y 2 y 4 serán iguales si c=d,                                                    
        {theta_ck, pg4_c_sols[1][0], pg4_c_sols[1][1]},    // y las soluciones 1 y 2 y 3 y 4 serán iguales si los 
        {theta_dk, pg4_d_sols[0][0], pg4_d_sols[0][1]},    // puntos intermedios de pg4 son iguales. Si c=d y los puntos intermedios
        {theta_dk, pg4_d_sols[1][0], pg4_d_sols[1][1]}     // de pg4 también, las cuatro soluciones serán iguales
    };

    return true;
}

// -----------------------------------------------------------------------------