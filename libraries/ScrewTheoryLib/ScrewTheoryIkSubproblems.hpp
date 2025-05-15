// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

#ifndef __SCREW_THEORY_IK_SUBPROBLEMS_HPP__
#define __SCREW_THEORY_IK_SUBPROBLEMS_HPP__

#include <kdl/frames.hpp>

#include "ScrewTheoryIkProblem.hpp"
#include "MatrixExponential.hpp"

namespace roboticslab
{

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief First Paden-Kahan subproblem
 *
 * Single solution, single revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot p = k @f$
 * (rotation screw applied to a point).
 */
class PadenKahanOne : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p Characteristic point.
     */
    PadenKahanOne(const MatrixExponential & exp, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PK1"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Second Paden-Kahan subproblem
 *
 * Dual solution, double revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive crossing rotation screws to a point).
 */
class PadenKahanTwo : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     * @param r Point of intersection between both screw axes.
     */
    PadenKahanTwo(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p, const KDL::Vector & r);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PK2"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, r, axesCross;
    const KDL::Rotation axisPow1, axisPow2;
    const double axesDot;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Third Paden-Kahan subproblem
 *
 * Dual solution, single revolute joint geometric IK subproblem given by
 * @f$ \left \| e\,^{\hat{\xi}\,{\theta}} \cdot p - k \right \| = \delta @f$
 * (rotation screw for moving @f$ p @f$ to a distance @f$ \delta @f$ from @f$ k @f$).
 */
class PadenKahanThree : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p First characteristic point.
     * @param k Second characteristic point.
     */
    PadenKahanThree(const MatrixExponential & exp, const KDL::Vector & p, const KDL::Vector & k);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PK3"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p, k;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief First Pardos-Gotor subproblem
 *
 * Single solution, single prismatic joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi}\,{\theta}} \cdot p = k @f$
 * (translation screw applied to a point, see @cite pardosgotor2018str
 * @cite pardosgotor2022str).
 */
class PardosGotorOne : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p Characteristic point.
     */
    PardosGotorOne(const MatrixExponential & exp, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PG1"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Second Pardos-Gotor subproblem
 *
 * Single solution, double prismatic joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive translation screws to a point, see @cite pardosgotor2018str
 * @cite pardosgotor2022str).
 */
class PardosGotorTwo : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     */
    PardosGotorTwo(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 1; }

    const char * describe() const override
    { return "PG2"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, crossPr2;
    const double crossPr2Norm;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Third Pardos-Gotor subproblem
 *
 * Dual solution, single prismatic joint geometric IK subproblem given by
 * @f$ \left \| e\,^{\hat{\xi}\,{\theta}} \cdot p - k \right \| = \delta @f$
 * (translation screw for moving @f$ p @f$ to a distance @f$ \delta @f$ from @f$ k @f$,
 * see @cite pardosgotor2018str @cite pardosgotor2022str).
 */
class PardosGotorThree : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp POE term.
     * @param p First characteristic point.
     * @param k Second characteristic point.
     */
    PardosGotorThree(const MatrixExponential & exp, const KDL::Vector & p, const KDL::Vector & k);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG3"; }

private:
    const MatrixExponential exp;
    const KDL::Vector p, k;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Fourth Pardos-Gotor subproblem
 *
 * Dual solution, double revolute joint geometric IK subproblem given by
 * @f$ e\,^{\hat{\xi_1}\,{\theta_1}} \cdot e\,^{\hat{\xi_2}\,{\theta_2}} \cdot p = k @f$
 * (consecutive parallel rotation screws applied to a point,
 * see @cite pardosgotor2018str @cite pardosgotor2022str).
 */
class PardosGotorFour : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term.
     * @param exp2 Second POE term.
     * @param p Characteristic point.
     */
    PardosGotorFour(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG4"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, n;
    const KDL::Rotation axisPow;
};

/**
 * @ingroup ScrewTheoryLib
 *
 * @brief Sixth Pardos-Gotor subproblem
 *
 * Dual solution, double revolute joint geometric IK subproblem defined as:
 * @f$ e^{\hat{\xi}_1 \theta_1} e^{\hat{\xi}_2 \theta_2} p = k @f$
 * (two consecutive skew rotation screws applied to a point).
 * This is a generic subproblem (PG6), generalizing PK2 and PG4.
 *
 * It handles the case of two skew axes, computes intermediate positions
 * using a translational screw axis (intersection of planes), and solves
 * two PK1 subproblems to get @f$ \theta_1 @f$ and @f$ \theta_2 @f$.
 *
 * See @cite pardosgotor2018str @cite pardosgotor2022str.
 */
class PardosGotorSix : public ScrewTheoryIkSubproblem
{
public:
    using ScrewTheoryIkSubproblem::solve;

    /**
     * @brief Constructor
     *
     * @param exp1 First POE term (rotation screw).
     * @param exp2 Second POE term (rotation screw).
     * @param p Initial point.
     * @param k Target point.
     */
    PardosGotorSix(const MatrixExponential & exp1, const MatrixExponential & exp2, const KDL::Vector & p);

    bool solve(const KDL::Frame & rhs, const KDL::Frame & pointTransform, const JointConfig & reference, Solutions & solutions) const override;

    int solutions() const override
    { return 2; }

    const char * describe() const override
    { return "PG6"; }

private:
    const MatrixExponential exp1, exp2;
    const KDL::Vector p, k;
    const KDL::Rotation axisPow1, axisPow2;

    static KDL::Vector computeR3(const KDL::Vector & o1, const KDL::Vector & o2,
                                 const KDL::Vector & omega1, const KDL::Vector & omega2);

    static std::array<double, 2> computeTheta3(const KDL::Vector & v3, const KDL::Vector & o,
                                               const KDL::Vector & r3, double targetNorm);

    static std::array<double, 2> solvePK1(const KDL::Vector & omega, const KDL::Vector & u,
                                          const KDL::Vector & m, const KDL::Vector & n);
};


} // namespace roboticslab

#endif // __SCREW_THEORY_IK_SUBPROBLEMS_HPP__
