#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memory"
#include "module_base/mathzone.h"
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include <vector>
#include <valarray>
#ifdef __MPI
#include "mpi.h"
#endif
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
#endif
Magnetism::Magnetism()
{
	this->tot_magnetization = 0.0;
	this->abs_magnetization = 0.0;
	this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
	delete[] this->start_magnetization;
}

/************************************************
 *  unit test of class UnitCell
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadNaoFlzTest
 *     - read_nao_flz(): read numerical orbital file and store in atom member variables rgrid, flz and n_rgrid
 * - NOTE: This function is only a templorary choice. Once module_nao is decoupled from old class ORB, this function will be removed.
 *         Currently the function to test is for reading-in detailed numerical orbital information for so-called lcao_in_pw calculation
 *         but now which is actually, the initialization by numerical orbitals on psi, in planewave representation
 * - To remove:
 *   - remove definition and declaration in read_atoms.cpp and unitcell.h
 *   - remove present test in CMakeLists.txt
 *   - delete this file
 */

//mock function
#ifdef __LCAO
void LCAO_Orbitals::bcast_files(
	const int &ntype_in,
	const int &my_rank)
{
	return;
}
#endif

class UcellTest : public ::testing::Test
{
protected:
    /* Here mock one unitcell object */
	UcellTestPrepare utp = UcellTestLib["flz-Read"];
	std::unique_ptr<UnitCell> ucell;
	std::ofstream ofs;
	std::string pp_dir;
	std::string output;
	void SetUp()
	{
		ofs.open("running.log");
		GlobalV::relax_new = utp.relax_new;
		GlobalV::global_out_dir = "./";
		ucell = utp.SetUcellInfo();
		GlobalV::LSPINORB = false;
		pp_dir = "./support/";
		GlobalV::PSEUDORCUT = 15.0;
		GlobalV::DFT_FUNCTIONAL = "default";
		GlobalV::test_unitcell = 1;
		GlobalV::test_pseudo_cell = 1;
		GlobalV::NSPIN = 1;
		GlobalV::BASIS_TYPE = "pw";

        ucell->atoms[0].l_nchi = new int[3];
        ucell->atoms[0].nwl = 2;
        ucell->atoms[0].l_nchi[0] = 2;
        ucell->atoms[0].l_nchi[1] = 2;
        ucell->atoms[0].l_nchi[2] = 1;

        ucell->atoms[1].l_nchi = new int[2];
        ucell->atoms[1].nwl = 1;
        ucell->atoms[0].l_nchi[0] = 2;
        ucell->atoms[0].l_nchi[1] = 1;
	}
	void TearDown()
	{
		ofs.close();
	}
};

using UcellDeathTest = UcellTest;

TEST_F(UcellDeathTest, ReadNaoFlzTest)
{
    std::string fn = "C_gga_8au_100Ry_2s2p1d.orb";
    ucell->read_nao_flz(0, fn, ofs, &(ucell->atoms[0]));
    EXPECT_EQ(ucell->atoms[0].n_rgrid, 801);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[0][0], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[0][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[0][800], 8);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[1][0], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[1][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].rgrid[1][800], 8);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[0][0], 5.368426038998e-01);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[0][4], 5.386156949024e-01);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[0][800], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[1][0], -6.134205291735e-02);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[1][4], -5.862821438615e-02);
    EXPECT_DOUBLE_EQ(ucell->atoms[0].flz[1][800], 0);
    fn = "H_gga_8au_100Ry_2s1p.orb";
    ucell->read_nao_flz(1, fn, ofs, &(ucell->atoms[1]));
    EXPECT_EQ(ucell->atoms[1].n_rgrid, 801);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[0][0], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[0][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[0][800], 8);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[1][0], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[1][1], 0.01);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].rgrid[1][800], 8);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[0][0], 1.857864679053e+00);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[0][4], 1.853731330136e+00);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[0][800], 0);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[1][0], 2.519718830629e+00);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[1][4], 2.511689716984e+00);
    EXPECT_DOUBLE_EQ(ucell->atoms[1].flz[1][800], 0);
}