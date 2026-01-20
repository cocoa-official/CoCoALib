//   Copyright (c)  2022-2023  John Abbott and Anna M. Bigatti
//   Code developed from original CoCoA-3 code by Fabio Rossi (1999)
//   Major contributions from Nico Mexis (2022-2023)

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/SparsePolyOps-cyclotomic.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingElemOps-CoprimeFactorBasis.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-eval.H"
#include "CoCoA/verbose.H"

#include "CoCoA/time.H"

//#include <functional>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {


    //------------------------------------------------------------------
    // Coefficient height bounds for cyclotomic polynomials
    // Various sorts of bound

    // // The bounds below were taken from a poster by Arnold+Monagan (CECM 2009)
    // // Result is a prime modulus greater than 2*CoeffHeight(cyclotomic(n)).
    // // Here index is the "index" of the cyclotomic poly.
    // long long cyclotomic_modulus(long index)
    // {
    //   if (index < 40755L) return 127LL;
    //   if (index < 327845L) return 2371LL;
    //   if (index < 707455L) return 62039LL;
    //   if (index < 1181895L) return 119633LL;
    //   if (index < 10163195L) return 150000001LL;
    //   if (index < 40324935L) return 5398416817471LL;
    //   CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    // }


//     /**
//      * Returns the highest coefficient a cyclotomic polynomial
//      * of squarefree index and given degree could possibly have.
//      * The comment after the line corresponds to the cyclotomic
//      * polynomial with the smallest degree that breaks the
//      * respective case.
//      */
//     long CycloCoeffHeightBound_sqfr(long d)
//     {
//       CoCoA_ASSERT(d > 1);
//       if (d == 1) return 1;
//       if (IsOdd(d)) return 0;
//       if (d%8 != 0) return 1; // at most 2 distinct primes [!!CHECK!!  well-known]
// //      if (d%16 == 8) return (2.0/3.0)*cbrt(d);  // at most 3 distinct primes "corrected sister Beiter conjecture"  SEE  arxiv:0910.2770.  Needs fixing?

//       // Values below here have been generated using a
//       // brute-force approach and are definitely correct.
//       if (d < 48) return 1; // Phi(105)
//       if (d < 240) return 2; // Phi(385)
//       if (d < 576) return 3; // Phi(1365)
//       if (d < 768) return 4; // Phi(1785)
//       if (d < 1280) return 5; // Phi(2805)
//       if (d < 1440) return 6; // Phi(3135)
//       if (d < 3840) return 7; // Phi(6545)
//       if (d < 5760) return 9; // Phi(15015)
//       if (d < 8640) return 23; // Phi(21945)
//       if (d < 10368) return 25; // Phi(25935)
//       if (d < 10560) return 27; // Phi(26565)
//       if (d < 17280) return 59; // Phi(40755)
//       if (d < 50688) return 359; // Phi(106743)
//       if (d < 82944) return 397; // Phi(171717)
//       if (d < 92160) return 434; // Phi(255255)
//       // Values below here have been generated using data
//       // by Arnold and Monagan, see here:
//       // http://wayback.cecm.sfu.ca/~ada26/cyclotomic/
//       // TODO: Data below is correct if and only if all heights are correct
//       if (d < 103680) return 532; // Phi(285285)
//       if (d < 126720) return 1182; // Phi(345345)
//       if (d < 138240) return 1311; // Phi(373065)
//       if (d < 193536) return 5477; // Phi(327845)
//       if (d < 387072) return 31010; // Phi(983535)
//       if (d < 483840) return 59518; // Phi(1181895)
//       if (d < 725760) return 14102773; // Phi(1752465)
//       if (d < 1824768) return 14703509; // Phi(3949491)
//       if (d < 3732480) return 56938657; // Phi(10555545)
//       if (d < 4147200) return 88835350; // Phi(11565015)
//       if (d < 4354560) return 197756850; // Phi(12267255)
//       if (d < 5806080) return 310102051; // Phi(10163195)
//       if (d < 7741440) return 1376877780831; // Phi(13441645)
//       if (d < 8709120) return 1475674234751; // Phi(15069565)
//       if (d < 11612160) return 1666495909761; // Phi(30489585)
//       // // TODO: Data below is most definitely incomplete
//       //   if (d < 15482880) return 2201904353336; // Phi(40324935)
//       //   if (d < 17418240) return 2699208408726; // Phi(43730115)
//       //   if (d < 104509440) return 862550638890874931; // Phi(306110805)
//       //   if (d < 231469056) return 4722828832054556497; // Phi(497111433)
//       //   if (d < 237828096) return 8171111062118177960; // Phi(516742863)
//       //   if (d < 240869376) return 8768227953282038629; // Phi(522080013)
//       //   if (d < 268240896) return 9038988754691465073; // Phi(582923523)
//       //   if (d < 321159168) return 9118090447189969651; // Phi(693722757)
//       //   if (d < 1072341504) return 9164566312887510757; // Phi(2583303555)
//       CoCoA_THROW_ERROR1(ERR::ArgTooBig);
//     }


    /**
     * Returns upper bound for coefficient height of a cyclotomic polynomial
     * of squarefree index <= k.
     * The comment after each line is the index of the first cyclotomic
     * polynomial whose coeff height is greater than all previous cyclos.
     */
    long CycloCoeffHeightBound_index(long k)
    {
      CoCoA_ASSERT(k > 0);
      CoCoA_ASSERT(MoebiusFn(k) != 0);
      if (IsEven(k))  k /= 2;
      if (k == 1)  return 1;
      // Values below here have been generated using a
      // brute-force approach and are definitely correct.
      if (k < 105)  return 1;
      if (k < 385)  return 2;
      if (k < 1365)  return 3;
      if (k < 1785)  return 4;
      if (k < 2805)  return 5;
      if (k < 3135)  return 6;
      if (k < 6545)  return 7;
      if (k < 10465)  return 9;
      if (k < 11305)  return 14;
      if (k < 17255)  return 23;
      if (k < 20615)  return 25;
      if (k < 26565)  return 27;
      if (k < 40755)  return 59;
      if (k < 106743)  return 359;
      if (k < 171717)  return 397;
      if (k < 255255)  return 434;
      // Values below here have been generated using data by Arnold and Monagan
      // see here:   http://wayback.cecm.sfu.ca/~ada26/cyclotomic/
      // TODO: Data below is correct if and only if all heights are correct
      if (k < 279565)  return 532;
      if (k < 327845)  return 1182;
      if (k < 707455)  return 31010;
      if (k < 886445)  return 35111;
      if (k < 983535)  return 44125;
      if (k < 1181895)  return 59518;
      if (k < 1752465)  return 14102773;
      if (k < 3949491)  return 14703509;
      if (k < 8070699)  return 56938657;
      if (k < 10163195)  return 74989473; // 32-bit limit!
#ifndef CoCoA_32BIT_LONG
      if (k < 13441645)  return 1376877780831;
      if (k < 15069565)  return 1475674234751;
      if (k < 30489585)  return 1666495909761;
      if (k < 40324935)  return 2201904353336;
      if (k < 43730115)  return 2699208408726;
      if (k < 306110805)  return 862550638890874931;
      if (k < 497111433)  return 4722828832054556497;
      if (k < 516742863)  return 8171111062118177960;
      if (k < 522080013)  return 8768227953282038629;
      if (k < 582923523)  return 9038988754691465073;
      if (k < 693722757)  return 9118090447189969651;
      if (k < 2583303555)  return 9164566312887510757; // 64-bit limit!
#endif
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    }


    // Returns const ref to a table of height bounds for first few coeffs of cyclo polys:
    // coeff of x^k has height at most tbl[k]
    const std::vector<int>& CyclotomicCoeffHeightTable()
    {
      // Table was computed using our own (inefficient) program.
      // Currently entries up to deg 180
      static vector<int> tbl{1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
        4, 3, 3, 3, 3, 4, 4, 5, 4, 4, 4, 5, 5, 6, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 9, 9, 7, 8, 8, 10,
        13, 12, 10, 12, 9, 11, 15, 13, 13, 14, 15, 13, 16, 15, 15, 14, 16, 24, 17, 21, 21, 16, 22,
        28, 26, 23, 28, 26, 25, 35, 34, 33, 28, 34, 36, 37, 49, 43, 33, 44, 48, 49, 55, 53, 53,
        48, 60, 70, 66, 65, 70, 65, 68, 91, 86, 78, 87, 86, 86, 109, 110, 98, 104, 108, 116, 124,
        136, 136, 118, 132, 153, 147, 162, 174, 150, 156, 187, 191, 196, 201, 194, 198, 213, 237,
        248, 229, 243, 254, 251, 294, 301, 291, 301, 316, 330, 337, 368, 375, 384, 393, 425, 434,
        444, 476, 495, 501, 534, 554, 575, 585, 633, 648, 659, 697, 727, 762, 770, 826, 839, 870,
        927, 939, 980, 1023, 1058, 1106, 1135, 1195};
      // The big table below was taken from OEIS  (A138474)  --- !!! REQUIRES MORE THAN 32-BITS !!!
//      static vector<int> tbl{1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 4, 5, 4, 4, 4, 5, 5, 6, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 9, 9, 7, 8, 8, 10, 13, 12, 10, 12, 9, 11, 15, 13, 13, 14, 15, 13, 16, 15, 15, 14, 16, 24, 17, 21, 21, 16, 22, 28, 26, 23, 28, 26, 25, 35, 34, 33, 28, 34, 36, 37, 49, 43, 33, 44, 48, 49, 55, 53, 53, 48, 60, 70, 66, 65, 70, 65, 68, 91, 86, 78, 87, 86, 86, 109, 110, 98, 104, 108, 116, 124, 136, 136, 118, 132, 153, 147, 162, 174, 150, 156, 187, 191, 196, 201, 194, 198, 213, 237, 248, 229, 243, 254, 251, 294, 301, 291, 301, 316, 330, 337, 368, 375, 384, 393, 425, 434, 444, 476, 495, 501, 534, 554, 575, 585, 633, 648, 659, 697, 727, 762, 770, 826, 839, 870, 927, 939, 980, 1023, 1058, 1106, 1135, 1195, 1222, 1270, 1345, 1364, 1420, 1475, 1551, 1576, 1648, 1723, 1765, 1833, 1917, 1975, 2025, 2128, 2205, 2277, 2354, 2449, 2533, 2609, 2716, 2822, 2887, 3033, 3118, 3214, 3352, 3463, 3577, 3684, 3852, 3959, 4091, 4249, 4403, 4515, 4695, 4879, 5013, 5166, 5384, 5557, 5711, 5931, 6140, 6313, 6547, 6783, 6968, 7215, 7444, 7729, 7927, 8230, 8499, 8735, 9068, 9350, 9640, 9940, 10311, 10611, 10917, 11333, 11691, 12006, 12442, 12831, 13240, 13619, 14106, 14539, 14933, 15478, 15967, 16423, 16934, 17520, 18004, 18586, 19166, 19754, 20342, 20998, 21669, 22265, 22980, 23711, 24386, 25136, 25935, 26681, 27474, 28348, 29239, 29988, 30974, 31939, 32821, 33802, 34849, 35870, 36885, 38040, 39161, 40280, 41450, 42768, 43926, 45192, 46604, 47899, 49276, 50765, 52212, 53706, 55264, 56912, 58451, 60136, 61982, 63636, 65476, 67336, 69337, 71177, 73261, 75421, 77388, 79611, 81943, 84198, 86479, 88977, 91506, 93969, 96621, 99375, 102015, 104883, 107807, 110812, 113731, 116995, 120213, 123396, 126862, 130352, 133832, 137440, 141346, 145032, 148957, 153094, 157208, 161350, 165736, 170297, 174767, 179350, 184367, 189151, 194125, 199442, 204712, 210066, 215629, 221503, 227152, 233159, 239373, 245635, 251993, 258727, 265467, 272263, 279482, 286802, 294149, 301692, 309723, 317641, 325697, 334313, 342879, 351490, 360736, 369949, 379351, 388924, 399069, 409187, 419345, 430271, 441119, 452211, 463627, 475516, 487284, 499526, 512291, 524978, 538047, 551691, 565467, 579416, 593915, 608802, 623715, 639190, 655254, 671227, 687736, 704833, 722346, 739702, 758069, 776794, 795529, 815042, 835084, 855346, 875912, 897527, 919313, 941351, 964230, 987691, 1011324, 1035647, 1060738, 1086214, 1111981, 1139041, 1166215, 1193792, 1222525, 1251751, 1281321, 1311707, 1343183, 1374747, 1407242, 1440661, 1474827, 1509149, 1544916, 1581526, 1618225, 1656179, 1695286, 1734764, 1775063, 1816717, 1859164, 1902034, 1946402, 1991886, 2037459, 2084865, 2133122, 2182387, 2232360, 2284047, 2336720, 2389897, 2445029, 2501160, 2557987, 2616386, 2676613, 2737215, 2799275, 2863350, 2928336, 2994302, 3062371, 3131806, 3202293, 3274266, 3348628, 3423647, 3500224, 3579361, 3659521, 3741172, 3824800, 3910779, 3997378, 4086497, 4177827, 4270325, 4364969, 4462138, 4560960, 4661336, 4764549, 4870016, 4976809, 5086277, 5198639, 5312379, 5428620, 5547998, 5669559, 5792537, 5919572, 6048904, 6180016, 6314331, 6451990, 6591777, 6734099, 6880555, 7029151, 7180596, 7335548, 7494178, 7654735, 7819038, 7987714, 8158545, 8332878, 8511885, 8693506, 8878694, 9068287, 9261752, 9458090, 9658820, 9864874, 10073174, 10286286, 10504265, 10726457, 10951921, 11183197, 11419241, 11658441, 11903358, 12153908, 12408233, 12667296, 12933110, 13203358, 13478069, 13759405, 14046293, 14337627, 14635630, 14939674, 15249367, 15564087, 15887182, 16215398, 16548944, 16890901, 17238932, 17592887, 17954180, 18323771, 18698613, 19081147, 19472467, 19870489, 20275020, 20689118, 21111202, 21539918, 21977359, 22424867, 22879065, 23341933, 23815498, 24296999, 24786837, 25287266, 25798069, 26315933, 26845537, 27385533, 27934826, 28494300, 29065902, 29647613, 30238885, 30843564, 31459281, 32084985, 32723423, 33375450, 34037643, 34712117, 35401525, 36102674, 36815222, 37544138, 38285512, 39039548, 39808414, 40593184, 41390912, 42202360, 43032131, 43875571, 44733722, 45609136, 46501548, 47407781, 48332455, 49275586, 50233637, 51209627, 52205803, 53218943, 54249394, 55300736, 56371639, 57459683, 58569272, 59700850, 60849879, 62021011, 63214983, 64430526, 65664955, 66925458, 68209042, 69512466, 70841932, 72196640, 73573572, 74975144, 76405173, 77859037, 79337725, 80845117, 82380574, 83940465, 85529872, 87149621, 88796266, 90471567, 92180778, 93917689, 95684436, 97485841, 99318840, 101182135, 103079538, 105013981, 106978275, 108978894, 111017058, 113090028, 115197382, 117345431, 119531909, 121752840, 124015115, 126320150, 128661962, 131044841, 133473237, 135942371, 138452614, 141010013, 143613043, 146256668, 148950464, 151691862, 154479303, 157314064, 160202182, 163139008, 166123741, 169164825, 172258011, 175401528, 178601623, 181860409, 185170345, 188538373, 191968832, 195455707, 198999885, 202609804, 206281142, 210012687, 213809344, 217675255, 221602541, 225596648, 229665200, 233799176, 238002237, 242279736, 246632470, 251053432, 255553048, 260131885, 264784053, 269515832, 274332266, 279227482, 284202775, 289267184, 294417097, 299649689, 304973052, 310389923, 315892862, 321488906, 327183894, 332972515, 338852949, 344840325, 350926331, 357109231, 363399451, 369797500, 376298953, 382906955, 389632406, 396465824, 403410376, 410475327, 417658877, 424955348, 432376454, 439925092, 447593278, 455387678, 463318079, 471375012, 479562708, 487891027, 496356256, 504955762, 513700002, 522594052, 531625490, 540807781, 550145213, 559634207, 569273280, 579077070, 589042114, 599162068, 609452787, 619915935, 630542474, 641341995, 652325455, 663483188, 674818298, 686343654, 698057289, 709953136, 722047614, 734340136, 746827701, 759514418, 772415624, 785520383, 798830568, 812364948, 826115567, 840081563, 854275788, 868704223, 883354993, 898243016, 913376660, 928748327, 944360399, 960231307, 976356018, 992729896, 1009369098, 1026280778, 1043452831, 1060898413, 1078630767, 1096639281, 1114929271, 1133516627, 1152402460, 1171574401, 1191059041, 1210855439, 1230957553, 1251377843, 1272127778, 1293201452, 1314600753, 1336347039, 1358435589, 1380861992, 1403646456, 1426796564, 1450298798, 1474170462, 1498426087, 1523056079, 1548063636, 1573475540, 1599281087, 1625481551, 1652096091, 1679132681, 1706581303, 1734453900, 1762774705, 1791526948, 1820721145, 1850376366, 1880494106, 1911066421, 1942119404, 1973660510, 2005677481, 2038189422, 2071215942, 2104743966, 2138782604, 2173356586, 2208462633, 2244098661, 2280288462, 2317043119, 2354347491, 2392227481, 2430697903, 2469753457, 2509394662, 2549657544, 2591776280, 2638859718, 2686767003, 2735512622, 2785109612, 2835572150, 2886916559, 2939154049, 2992304101, 3046379312, 3101394489, 3157368382, 3214314131, 3272249082, 3331191327, 3391154353, 3452157970, 3514218830, 3577353722, 3641582150, 3706921226, 3773388432, 3841005078, 3909789193, 3979758493, 4050935817, 4123338326, 4196988118, 4271905516, 4348110945, 4425625930, 4504473234, 4584672936, 4666249343, 4749224187, 4833620167, 4919461259, 5006772958, 5095575512, 5185897557, 5277761714, 5371193275, 5466221231, 5562866803, 5661159913, 5761127462, 5862795576, 5966192555, 6071348757, 6178288470, 6287045682, 6397647790, 6510124309, 6624508040, 6740827957, 6859116008, 6979407790, 7101729658, 7226119165, 7352609762, 7481233156, 7612027191, 7745024551, 7880261728, 8017775831, 8157603970, 8299781531, 8444349539, 8591344032, 8740805696, 8892775513, 9047291277, 9204395314, 9364131129, 9526538189, 9691662422, 9859545604, 10030232379, 10203768784, 10380201948, 10559573835, 10741936359, 10927335890, 11115819575, 11307441081, 11502245790, 11700287269, 11901618476, 12106288801, 12314354530, 12525870093, 12740885971, 12959465157, 13181659717, 13407527625, 13637129315, 13870523446, 14107768700, 14348931413, 14594066787, 14843242863, 15096524137, 15353971519, 15615656396, 15881642451, 16151997998, 16426794165, 16706100821, 16989987410, 17278529374, 17571796964, 17869867308, 18172816465, 18480719655, 18793654490, 19111704711, 19434945061, 19763461388, 20097335722, 20436650194, 20781493763, 21131952207, 21488111319, 21850063877, 22217898864, 22591707549, 22971588230, 23357630341, 23749931570, 24148594093, 24553711592, 24965388342, 25383727867, 25808827250, 26240802067, 26679752478, 27125788226, 27579021974, 28039563714, 28507526546, 28983030370, 29466185791, 29957115645, 30455941475, 30962781835, 31477765659, 32001016258, 32532659582, 33072831591, 33621660953, 34179278767, 34745826645, 35321437086, 35906253977, 36500417687, 37104072809, 37717362928, 38340442814, 38973456404, 39616560551, 40269910405, 40933659594, 41607973273, 42293011810, 42988936567, 43695919518, 44414126897, 45143730656, 45884910118, 46637835225, 47402688872, 48179657390, 48968919455, 49770665819, 50585090475, 51412377408, 52252734311, 53106353128, 53973435992, 54854192944, 55748826468};
      return tbl;
    }




    // --------------------------------------------
    // Compute cyclo prefix using the Moebius inverted product
    // We implement ad hoc truncated power-series (as vector<long>)

    // !!! SUGGESTION !!!
    // Make this into a template fn: replace "unsigned long" by T;
    // change only type decl of cyclo, and return type.
    // Then we can instantiate with unsigned int types or BigInt!

    // Dmax > 0 means Compute the first Dmax+1 coeffs of the cyclotomic poly (in a std::vector)
    // Dmax == 0 implies use half the degree;
    // Dmax < 0 means compute the abs(Dmax)+1 prefix of RECIPROCAL of the requested cyclo
    std::vector<unsigned long> CycloPrefix(const std::vector<long>& plist,  long Dmax=0)
    {
      CoCoA_ASSERT(!plist.empty());  // entries are DISTINCT PRIMES -- NOT CHECKED!!!

      vector<long> even;
      vector<long> odd;
      even.push_back(1);
      long phi = 1;
      bool Prime2Present = false;
      // This loop computes all products of even/odd subsets of plist; prime 2 is "ignored"
      for (int i=0; i < len(plist); ++i)
      {
        const long p = plist[i];
        if (Dmax != 0 && p > std::abs(Dmax))  continue;
        if (p == 2)  { Prime2Present = true; continue; }
        phi *= (p-1);
        const long UPB = std::abs(Dmax)/p; // we can ignore multiples of p by anything larger than UPB
        const int len_even = len(even);
        const int len_odd = len(odd);
        for (int j=0; j < len_even; ++j)
          if (Dmax==0 || even[j] <= UPB)
            odd.push_back(p*even[j]);
        for (int j=0; j < len_odd; ++j)
          if (Dmax==0 || odd[j] <= UPB)
            even.push_back(p*odd[j]);
      }
      const bool SwapOver = (Dmax < 0) ^ Prime2Present ^ IsOdd(len(plist));
      if ( SwapOver )
        swap(odd, even);
      // Now "even" corr to numerator; "odd" corr to denominator
      // NOTE: all entries of even & odd are <= abs(Dmax)
      sort(even.begin(), even.end());
      const long UPB = (Dmax==0) ? phi/2 : std::abs(Dmax); //??? or maybe std::min(std::abs(Dmax),phi/2)
      vector<unsigned long> cyclo(UPB+1); // unsigned important to avoid UB upon overflow
      cyclo[0] = 1;
      long d = 0; // current deg, but at most UPB
      // This loop computes numerator product:
      for (const long i: even)
      {
        // mult by (1+x^i)
        long top = d;
        if (top+i > UPB)
          top = UPB-i; // >= 0 by check above
        for (long j=top; j >= 0; --j)
          cyclo[j+i] += cyclo[j];
        d = std::min(d+i, UPB); /// d += i; if (d > UPB) d = UPB; // EQUIV:  d = min(d+i, UPB);
      }
      // Map x |-> -x
      for (long i=1; i <= UPB; i += 2)
        cyclo[i] = -cyclo[i];

      // This loop "divides" by denominator, one factor at a time.
      // Division is multiplication by truncated inverse power series:
      // e.g. 1/(1-x) = 1+x+x^2+x^3+...
      for (const long i: odd)
      {
        // div by (1-x^i)
        const long top = UPB - i; // >= 0 since all entries in odd are <= UPB
        for (long j=0; j <= top; ++j)
          cyclo[j+i] += cyclo[j];
      }
      if (Prime2Present)
      {
        // Map x |-> -x
        for (long i=1; i <= UPB; i += 2)
          cyclo[i] = -cyclo[i];
      }
      return cyclo;
    }




//     // Compute coeff vec of cyclotomic(n) mod modulus, where n = product(plist)
//     // plist is a list of odd primes in incr order.
//     // Let v be the result then v[k] is coeff of x^k in cyclotomic(n) mod modulus.
//     // Entries in v are symmetric remainders.
//     std::vector<long> cyclotomic_modp(long modulus, const vector<long>& plist)
//     {
//       // assume plist contains ODD primes in incr order
//       if (plist.empty()) return vector<long>(2,1);
//       long n=1; for (long p: plist) n *= p;
//       SmallFpImpl ModP(modulus);
//       DUPFp OnePoly(0,ModP);  AssignOne(OnePoly);
//       DUPFp tmp(n,ModP);

//       long p = plist[0];
// //    for (long k=0; k < p; ++k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
//       for (long k=p-1; k >= 0; --k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));

//       DUPFp phi = tmp;
//       long i=1;
//       while (i < len(plist))
//       {
//         CheckForInterrupt("cyclotomic_modp: gcd loop");
//         const long p = plist[i];
//         AssignZero(tmp);
// //      for (long k=0; k < p; ++k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
//         for (long k=p-1; k >= 0; --k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
//         phi = gcd(phi, tmp); // may take a long time
//         ++i;
//       }
//       const long d = deg(phi);
//       vector<long> ans(1+d);
//       for (long i=0; i <= d; ++i)
//         ans[i] = ModP.myExportSymm(phi.myCoeffs[i]);
//       return ans;
//     }



//     // Naive impl of CRT reconstruction -- probably good enough here
//     long CRT2(long r1, long m1, long r2, long m2)
//     {
//       long m1inv = InvMod(m1,m2);
//       long m2inv = InvMod(m2,m1);
//       long R = m2*SymmRemainder(r1*m2inv, m1) + m1*SymmRemainder(r2*m1inv, m2);
//       return SymmRemainder(R, m1*m2);
//     }



//   // Impl below is "curious": it computes the coeffs as machine "long int"
//   // then maps the result into the given ring.

//   // Developed from CoCoA code by Fabio Rossi (date 1999?)
//   RingElem cyclotomic_FabioRossi(long n, const RingElem& x)
//   {
//     if (n < 1)
//       CoCoA_THROW_ERROR1(ERR::ReqPositive);
//     if (!IsIndet(x))
//       CoCoA_THROW_ERROR1(ERR::ReqIndet);
//     if (n == 1)  return x-1;

//     int power2 = 0;
//     while (IsEven(n)) { ++power2; n /= 2; }
//     const factorization<long> facs = factor(n);
//     long radn = 1; for (long j: facs.myFactors())  { radn *= j; }
//     const long deg = EulerTotient(radn);
//     const vector<long> plist = facs.myFactors();
//     const long long M = cyclotomic_modulus(radn);
//     long p =  (M < 66000) ? M : NextPrime(1024);
//     long long CRT_modulus = 1;
//     vector<long> CoeffVec; // no need to reserve space
//     while (M > CRT_modulus)
//     {
//       auto CoeffVec_modp = cyclotomic_modp(p, plist);
//       if (CRT_modulus == 1)
//       { // here only for the 1st iter
//         CRT_modulus = p;
//         swap(CoeffVec, CoeffVec_modp); /*swap is assignment*/;
//         p = NextPrime(p);
//         continue;
//       }
//       // Here every iter after the first one
//       for (int i=0; i <= deg; ++i)
//       {
//         CoeffVec[i] = CRT2(CoeffVec[i], CRT_modulus, CoeffVec_modp[i],p);
//       }
//       CRT_modulus *= p;
//       p = NextPrime(p);
//     }

//     if (power2 > 0)
//     {
//       for (int i=1; i < deg; i += 2)
//         CoeffVec[i] = -CoeffVec[i];
//     }

//     // Convert result into a poly in ring of x
//     const long xpower = (power2 <= 1) ? (n/radn) : ((n/radn) << (power2-1));
//     const ring& P = owner(x);
//     const ring& R = CoeffRing(P);
//     const PPMonoidElem X = power(LPP(x), xpower);
//     RingElem ans(P);
//     for (int i=0; i <= deg; ++i)
//     {
//       if (CoeffVec[i] == 0) continue;
//       PushFront(ans, RingElem(R,CoeffVec[i]), power(X,i));
//     }
//     return ans;
//   }

  } // end of namespace anonymous




  // Calls CycloPrefix then builds the actual poly.
  // TODO make a version which uses CRT to allow larger indices?
  RingElem CyclotomicPoly(long n, ConstRefRingElem x)
  {
    if (n < 1)
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    long x_index;
    if (!IsIndet(x_index, x))
      CoCoA_THROW_ERROR1(ERR::ReqIndet);
    if (n == 1)  return x-1;
    const long OddPart = CoprimeFactor(n,2);
    if (OddPart == 1)  return power(x,n/2) + 1;
    const bool EvenIndex = IsEven(n);
    const long EvenPart = n/OddPart; // power of 2
    const factorization<long> facs = factor(OddPart);
    const long nfacs = len(facs.myFactors());
    long PhiOddPrimes = 1;
    n = OddPart;
    long OddRadical = 1;
    for (long p: facs.myFactors())
    { OddRadical *= p; PhiOddPrimes *= (p-1); } // neither prod can overflow!
    n /= OddRadical;
    // Safeguard against coeff overflow (since we use long int internally)
    if (nfacs > 3)
    {
      // Check only if nfacs > 3 because o/w coeffs have to be small
      if (sizeof(unsigned long) == 4 && OddRadical >= 10163195L)
        CoCoA_THROW_ERROR1(ERR::ArgTooBig);
      if (sizeof(unsigned long) == 8 && OddRadical >= 169828113L)
        CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    }

    const auto coeff = CycloPrefix(facs.myFactors());
    const long UPB = PhiOddPrimes/2; // exact division
    const long exponent = EvenIndex ? (n*EvenPart/2) : n;
    const PPMonoidElem X = power(LPP(x), exponent);
    { /*not used*/ const PPMonoidElem Xpwr = power(X, UPB); /* might ExpTooBig trigger error*/ if (deg(Xpwr)/UPB != deg(X))  CoCoA_THROW_ERROR1(ERR::ExpTooBig); }
    const vector<long> Xexpv = exponents(X);
    const int nvars = len(Xexpv);
    const SparsePolyRing P = owner(x);
    const ring& k = CoeffRing(P);
    RingElem coef(k);
    vector<long> expv(nvars);
    PPMonoidElem Xpwr = one(owner(X));
    RingElem ans(P);
    for (long i=0; i <= UPB; ++i)
    {
      if (coeff[i] == 0) continue;
      long C = ULong2Long(coeff[i]);
      if (EvenIndex && IsOdd(i))
        C = -C;
      coef = C; expv[x_index] = i*exponent;
//      coef = C; for(int j=0;j<nvars;++j)expv[j]=i*Xexpv[j];
      P->myPushFront(raw(ans), raw(coef), expv);
    }

    for (long i=UPB-1; i >= 0; --i)
    {
      if (coeff[i] == 0) continue;
      const PPMonoidElem t = power(X, 2*UPB-i);
      long C = ULong2Long(coeff[i]);
      if (EvenIndex && IsOdd(i))
        C = -C;
//      coef = C; for(int j=0;j<nvars;++j)expv[j]=(2*UPB-i)*Xexpv[j];
      coef = C; expv[x_index] = (2*UPB-i)*exponent;
      P->myPushFront(raw(ans), raw(coef), expv);
    }

    return ans;
  }




  //--------------------------------------------
  // Fns to check whether a poly is cyclo


  // This is a "naughty/ugly" function: it does several things at once
  //  (a) checks cyclo by matching against "prefixes" (& do only this if !DoFullCheck)
  //  (b) checks that f is univariate, palindromic, below height bound H
  //  (c) checks for special cases...
  // Input poly f
  // Output:
  //  (a) if input is "obviously not cyclo" return  0
  //  (b) special case: index is p^k or 2*p^k
  //  (c) special case: f(x) = g(x^r) for some r > 1
  unsigned long CyclotomicTest(ConstRefRingElem f, bool DoFullCheck)
  {
    constexpr unsigned long DefinitelyNotCyclo = 0;  // "impossible" index
    if (IsConstant(f))  return DefinitelyNotCyclo;
    if (!IsOne(LC(f)))  return DefinitelyNotCyclo;
    const long degf = deg(LPP(f)); // assume even, monic, not of form g(x^r), and coeff of x^(degf-1) is +-1
    SparsePolyIter it = BeginIter(f);
    ++it;  // skip LPP -- already handled it above
    if (IsEnded(it))  return DefinitelyNotCyclo; // cannot be cyclo with just a single summand!
    // Do special case deg(f) == 1:
    if (degf == 1)
    {
      if (!IsOne(PP(it)))  return DefinitelyNotCyclo;
      if (IsMinusOne(coeff(it)))  return 1;
      if (IsOne(coeff(it)))  return 2;
      return DefinitelyNotCyclo;
    }
    if (IsOdd(degf))  return DefinitelyNotCyclo;
    const PPMonoidElem x = radical(LPP(f));
    if (deg(x) != 1)
      CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    if (IsOne(PP(it)))
    {
      // Poly is of form x^k+nzconst
      if (IsOne(coeff(it)) && CoprimeFactor(degf,2) == 1)
        return 2*degf; // ???BUG??? overflow???
      else
        return DefinitelyNotCyclo;
    }

    const long r = degf - deg(PP(it));
    if (r <= 0) // cannot be univariate
      CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    if (degf%r != 0)  return DefinitelyNotCyclo; // not of form g(x^r)
    if (!IsOne(coeff(it)) && !IsMinusOne(coeff(it)))  return DefinitelyNotCyclo;
    const long mu = (IsOne(coeff(it)))? -1 : 1; // minus 2nd coeff
    const long radr = radical(r);
    const long degr = degf/r; // "reduced degree"
    // Obtain list of candidate indices to try: in 3 stages
    // InvTotient (sqfr preimages), with correct MoebiusFn value, & must be mult of r
    vector<long> cand = InvTotient(degr, InvTotientMode::SqFreePreimages);
    const auto WrongMu = [mu](long n){ return (MoebiusFn(n) != mu); };
    { auto it = std::remove_if(cand.begin(), cand.end(), WrongMu);  cand.erase(it, cand.end()); } // C++20 erase_if(...)
    const auto NotMultOfRadr = [radr](long n){return (n%radr != 0);};
//C++20    if (radr != 1)  erase_if(cand, [radr](long n){return (n%radr != 0);});
    if (radr != 1)
    { auto it = std::remove_if(cand.begin(), cand.end(), NotMultOfRadr);  cand.erase(it, cand.end()); }
    if (cand.empty())  return DefinitelyNotCyclo;
    const vector<int>& CycloCoeffHeightTbl = CyclotomicCoeffHeightTable();
    vector<long> CandHeightBound;
    long H=0;  // overall height bound -- 0 means "to be computed"
    vector<long> C(1+degr/2);  C[0] = 1;
    long prev = 2;
    long thresh = degr/2;  while (thresh >= 64)  { thresh = 1+thresh/4; }

    // This loop handles the "upper half"
    while (!IsEnded(it))
    {
      const long d = deg(PP(it));
      if (d < degf/2)  break;
      if (radical(PP(it)) != x)
        CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
      const long dr = (degf-d)/r;
      if (degf-d != r*dr)  return DefinitelyNotCyclo;
      long c;
      if (!IsConvertible(c, coeff(it)))  return DefinitelyNotCyclo;
      if (dr < len(CycloCoeffHeightTbl))
      { if (std::abs(c) > CycloCoeffHeightTbl[dr])  return DefinitelyNotCyclo; }
      else { if (H==0) { long kmax=0;for (long k: cand) kmax = std::max(kmax,k); H = CycloCoeffHeightBound_index(kmax); }
        if (std::abs(c) > H)  { return DefinitelyNotCyclo; } }
      C[dr] = c;
      ++it;
      if (IsEnded(it))  return DefinitelyNotCyclo;
      const long deg_next = deg(PP(it));
      if (deg_next%r != 0)  return DefinitelyNotCyclo;
      // OR MAYBE:   if (dr < thresh) continue;
      if ((degf-deg_next)/r >= thresh) //(dr >= thresh)
      {
        // Coeffs in C are correct up to index dr_next, and dr_next >= thresh:
        const long dr_next = (deg_next < degf/2) ? degf/(2*r) : (degf-deg_next)/r-1;
        vector<long> NewCand;
        for (long k: cand)
        {
          const factorization<long> facs = factor(k);
          const vector<unsigned long> prefix = CycloPrefix(facs.myFactors(), dr_next);
          bool eq = true;
          for (int i=prev; i <= dr_next; ++i)
            if (prefix[i] != static_cast<unsigned long>(C[i]))  { eq = false; break; }
          if (eq)  NewCand.push_back(k);
        }
        if (NewCand.empty())  return DefinitelyNotCyclo;
        if (len(cand) != len(NewCand)) { swap(cand, NewCand); H = 0/*force recomputation*/; }
        if (!DoFullCheck && len(cand) == 1 && dr_next > 31)  return r*cand[0]; // might be false positive!
        prev = dr_next+1;
        do {thresh *= 4;} while (thresh <= 2*dr_next);
        thresh = std::min(thresh, degr/2);
      }
    }
    if (IsEnded(it)) return DefinitelyNotCyclo;

    // Loop below checks that f is palindromic (actually a bit slow)
    long PredictedDegr = degr/2;
    while (!IsEnded(it))
    {
      while (C[--PredictedDegr] == 0); // SAFE!  PredictedDegr can never become negative since C[0] != 0
      const long d = deg(PP(it));
      if (d != r*PredictedDegr)  return DefinitelyNotCyclo;
      if (d > 0 && radical(PP(it)) != x)  // because not univariate
        CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
//      if (coeff(it) != C[PredictedDegr])  return DefinitelyNotCyclo;  // SLOWER THAN CODE BELOW
      long c;
      if (!IsConvertible(c,coeff(it)) || c != C[PredictedDegr])  return DefinitelyNotCyclo;
      ++it;
    }
    return r*cand[0];
  }




  //------------------------------------------------------------------
  // Impl of CyclotomicFactors (developed from Smyth+Beukers)

  namespace // anonymous
  {
    // Input is f(x^2); output is f(x)
    RingElem sqrtx(ConstRefRingElem f)
    {
      // Assume f is univariate;  non-zero???
      if (IsConstant(f)) return f;
      const ring& P = owner(f); // assume sparse poly ring
      PPMonoidElem x = radical(LPP(f));
      RingElem ans(P);
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        PushBack(ans, coeff(it), power(x,deg(PP(it))/2));  // assume exact division!
      }
      return ans;
    }


    void CONCAT_move(vector<RingElem>& v, vector<RingElem> v2)
    {
      // Taken from https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
      v.insert(
        v.end(),
        std::make_move_iterator(v2.begin()),
        std::make_move_iterator(v2.end())
               );
    }


    std::vector<RingElem> BeukersSmythOpG(ConstRefRingElem f, long Xindex)
    {
      VerboseLog VERBOSE("BeukersSmythOpG");
      vector<RingElem> ans;
      if (IsConstant(f))  return ans;
      // assume f is univariate (non-const?)
      const ring& P = owner(f);
      const long nvars = NumIndets(P);
      vector<RingElem> ImagesNegate(nvars, zero(P));
      ImagesNegate[Xindex] = -indet(P,Xindex);
      RingHom NegateX = PolyAlgebraHom(P,P,ImagesNegate);
      RingElem fneg = NegateX(f);
      VERBOSE(85) << "Computing even factor" << std::endl;
      RingElem f2 = gcd(f, fneg);
      if (IsConstant(f)) { VERBOSE(85) << "Trivial even factor" << std::endl; return ans; } // empty list
      vector<RingElem> ImagesSquare(nvars, zero(P));
      ImagesSquare[Xindex] = power(indet(P,Xindex),2);
      RingHom SquareX = PolyAlgebraHom(P,P, ImagesSquare);
      RingElem g1 = SquareX(f/f2);
      RingElem g2 = SquareX(f2);
      RingElem g3 = SquareX(fneg/f2);
      VERBOSE(85) << "Computing 3 gcds..." << std::endl;
      VERBOSE(89) << "GCD 1" << std::endl;
      VERBOSE(89) << "deg(f) = " << deg(LPP(f)) << std::endl;
      VERBOSE(89) << "deg(g1) = " << deg(LPP(g1)) << std::endl;
      VERBOSE(89) << "deg(g2) = " << deg(LPP(g2)) << std::endl;
      VERBOSE(89) << "deg(g3) = " << deg(LPP(g3)) << std::endl;
      RingElem h = gcd(f,g1);
      VERBOSE(89) << "gcd(f,g1) has deg " << deg(LPP(h)) << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(89) << "GCD 2" << std::endl;
      h = gcd(f,g2);
      VERBOSE(89) << "gcd(f,g2) has deg " << deg(LPP(h)) << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(89) << "GCD 3" << std::endl;
      h = gcd(f,g3);
      VERBOSE(89) << "gcd(f,g3) has deg " << deg(LPP(h)) << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(85) << "Finished" << std::endl;
      return ans;
    }

    std::vector<RingElem> BeukersSmythOpC(ConstRefRingElem f, long Xindex)
    {
      VerboseLog VERBOSE("BeukersSmythOpC");
      if (IsConstant(f))  return vector<RingElem>();
      VERBOSE(85) << "Applying OpG" << std::endl;
      vector<RingElem> L = BeukersSmythOpG(f, Xindex);
      if (!L.empty())
      {
        long d=0; for (const RingElem& g: L)  d += deg(g);
        if (d == deg(f))  { VERBOSE(85) << "Deg check short cut" << std:: endl; return L; }
      }
      // Copy-paste from above (redundant!!)
      const ring& P = owner(f);
      const long nvars = NumIndets(P);
      vector<RingElem> ImagesNegate(nvars, zero(P));
      ImagesNegate[Xindex] = -indet(P,Xindex);
      RingHom NegateX = PolyAlgebraHom(P,P,ImagesNegate);
      vector<RingElem> ImagesSquare(nvars, zero(P));
      ImagesSquare[Xindex] = power(indet(P,Xindex),2);
      RingHom SquareX = PolyAlgebraHom(P,P, ImagesSquare);

      VERBOSE(85) << "Checking even factor" << std::endl;
      const RingElem f2 = gcd(f, NegateX(f));
//      const RingElem g = sqrtx(f2);
      VERBOSE(85) << "Applying OpC to even factor" << std::endl;
      vector<RingElem> Lg = BeukersSmythOpC(sqrtx(f2), Xindex);
      for (int i=0; i < len(Lg); ++i)  Lg[i] = SquareX(Lg[i]);
      for (const RingElem& ff: L)  { CONCAT_move(Lg, BeukersSmythOpC(ff, Xindex)); }
      VERBOSE(85) << "Computing CoprimeFactorBasis" << std::endl;
      CoprimeFactorBasis_RingElem CFB;  CFB.myAddInfo(Lg);
      return FactorBase(CFB);
    }

  } // end of namespace anonymous



  // Returns a partial factorization of f:
  // each returned factor is a product of (distinct) cyclo factors of f;
  // the returned factors are pairwise coprime.
  factorization<RingElem> CyclotomicFactors_BeukersSmyth(ConstRefRingElem f)
  {
    const char* const FnName = "CyclotomicFactors_BeukersSmyth";
    VerboseLog VERBOSE(FnName);
    if (IsZero(f))
      CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (IsConstant(f))
      return factorization<RingElem>(f);
    const long x = UnivariateIndetIndex(f);
    if (x < 0)
      CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    factorization<RingElem> ans(one(owner(f)));
    VERBOSE(85) << "Extract palindromic factor" << std::endl;
    const RingElem pal = PalindromicFactor(f);
    if (IsConstant(pal))
      return ans; // empty factorization
    VERBOSE(85) << "Compute SquareFree factors" << std::endl;
    const auto facs = SqFreeFactor(pal);
    const vector<RingElem>& FacList = facs.myFactors();
    VERBOSE(85) << "Num square-free factors: " << len(FacList) << std::endl;
    for (int i=0; i < len(FacList); ++i)
    {
      const long m = facs.myMultiplicities()[i];
      VERBOSE(85) << "Doing factor " << i << " (mult=" << m << ", deg="<< deg(FacList[i]) << ")" << std::endl;
      for (const RingElem& g: BeukersSmythOpC(FacList[i],x))
        ans.myAppend(g, m);
    }
    return ans;
  }





// 2023-12-09 JAA: commented out this old version which took as input an "evaluator"
// //  CandidateCycloIndices is a list of integers >= 2 being indices of
// //  possible cyclotomic factors -- see also CycloIndicesUptoDeg (below).
// //  Deliberately skips indices 1 and 2
//     std::vector<long> FindCycloFactor(vector<long> CandidateCycloIndices, std::function<BigInt /**/(long n, long d)> EvalF)
//     {
//       // Assume f is non-constant with integer coeffs (ideally squarefree, content=1, and f palindromic)
//       VerboseLog VERBOSE("FindCycloFactor");

//       const ring ZZx = NewPolyRing(RingZZ(), symbols("x")); // used only to evaluate a cyclo poly
//       const RingElem& x = indet(ZZx,0);                    //


//       VERBOSE(80) <<"NUM INIT Candidates: " << len(CandidateCycloIndices) << std::endl;
//       int TargetNumEvalPts = 2;
//       long EvalPtNumer = 2;  long EvalPtDenom = 1;
//       int NumEvalPts = 0;
//       while (NumEvalPts < TargetNumEvalPts)
//       {
//         BigInt Valf = EvalF(EvalPtNumer, EvalPtDenom);
//         if (!IsZero(Valf))
//         {
//           ++NumEvalPts;
//           VERBOSE(80) <<"Chosen EvalPt="<<EvalPtNumer << "/" << EvalPtDenom <<std::endl;
//           BigInt Valf_reduced = CoprimeFactor(Valf, EvalPtNumer*(EvalPtNumer*EvalPtNumer-EvalPtDenom*EvalPtDenom)); // BUG?????  to avoid overflow need EvalPt^3 < MaxLong (or MaxULong)

//           vector<long> ReducedCandidateList;
//           for (long k:  CandidateCycloIndices)
//           {
//             if (k < 3)  continue; // skip 1 & 2, if present
//             if (EvalPtNumer == 2 && EvalPtDenom == 1 && k == 6 && Valf%3 == 0)
//             { ReducedCandidateList.push_back(k); continue; } // exception in Zsygmondi's Thm
//             /*const*/ BigInt g = power(EvalPtNumer, k) - power(EvalPtDenom,k);
//             g = gcd(g, Valf_reduced);
//             if (IsOne(g)) { /*clog << '*';*/ continue; }
//             Valf_reduced = CoprimeFactor(Valf_reduced, g);
// //        const BigInt ValCyclo = EvalAt(CyclotomicPoly(k,x),EvalPtNumer,EvalPtDenom);///ConvertTo<BigInt>(evaluate(CyclotomicPoly(k,x)));
//             if (/*EvalPtNumer != 2 ||*/ IsDivisible(Valf,EvalAt(CyclotomicPoly(k,x),EvalPtNumer,EvalPtDenom)))
//               ReducedCandidateList.push_back(k);
//             if (IsOne(Valf_reduced))  break;
//           }
//           if (len(ReducedCandidateList) < len(CandidateCycloIndices))
//           {
//             TargetNumEvalPts = NumEvalPts+3; // if no indices removed in 3 consecutive iters, we accept the list as "probably correct" (i.e. only very few false positives)
//             swap(CandidateCycloIndices, ReducedCandidateList); // really assignment
//             VERBOSE(80)<<"CandidateCycloIndices="<<CandidateCycloIndices<<std::endl;
//             if (CandidateCycloIndices.empty()) break;
//           }
//         }
//         // advance to next EvalPt
//         const long skip = RandomLong(1, 50); for (long i=0; i < skip; ++i)  GotoNextEvalPt(EvalPtNumer, EvalPtDenom);
//       }
//       VERBOSE(80)<<"Returning"<<std::endl;
//       return CandidateCycloIndices;
//     }

  namespace // anonymous -- auxiliaries for FindCycloFactor
  {

    // Update (n,d) to "next biggest rational" (in lex ordering)
    void  GotoNextEvalPt(long& n, long& d)
    {
      // assume n >= 2 and 1 <= d < n
      if (d >= n-2)  { ++n; d = 1; return; }  // trick to favour integers
      do { ++d; } while (gcd(n,d) != 1);
    }


    // struct used only inside FindCycloFactor (immediately below)
    struct SmallRat
    {
      SmallRat(unsigned short n, unsigned short d): num(n), den(d) {}
      unsigned short num;
      unsigned short den;
    };

  } // end of namespace anonymous


    //  CandidateCycloIndices is a list of integers >= 2 being indices of
    //  possible cyclotomic factors -- see also CycloIndicesUptoDeg (below).
    //  Deliberately skips indices 1 and 2
  // ASSUMES f is square-free
  std::vector<long> FindCycloFactor(vector<long> CandidateCycloIndices, RingElem f)
  {
    // Evaluation points such that low deg cyclos have a "large" prime factor

///    static SmallRat LastCheck[/*SizeLastCheck*/] = {{18,  17},  {21,  4}, {25,  7}, {25,  24}, {26,  5}, {26,  21}, {27,  13}, {31,  4}, {31,  10}, {31,  21}, {32,  3}, {33,  14}, {33,  20}, {34,  21}, {35,  8}, {35,  9}, {35,  18}, {35,  26}, {35,  31}, {35,  34}, {36,  13}, {36,  29}, {37,  2}, {37,  14}, {37,  15}, {37,  33}, {37,  35}, {38,  25}, {39,  10}, {39,  29}, {39,  32}, {39,  38}};
///    static SmallRat LastCheck[/*SizeLastCheck*/] = {{36, 29},  {39, 10},  {39, 35},  {40, 37},  {42, 5},  {45, 38},  {48, 43},  {49, 6},  {52, 45},  {55, 6},  {57, 8},  {61, 54},  {67, 52},  {70, 39},  {71, 6},  {72, 7},  {84, 19},  {88, 3},  {90, 1},  {96, 5},  {98, 3}};

//    static SmallRat LastCheck[/*SizeLastCheck*/] = {{111, 91},  {117, 98},  {117, 107},  {126, 121},  {129, 119},  {133, 18},  {133, 130},  {138, 5},  {144, 11},  {145, 142},  {147, 145},  {150, 7},  {159, 16},  {161, 159},  {165, 158},  {165, 164},  {169, 6},  {174, 161},  {175, 163},  {175, 172},  {183, 167},  {185, 18},  {187, 177},  {199, 24},  {200, 3}}; // cyclos 3,4,6 all have a large prime factor
    static SmallRat LastCheck[/*SizeLastCheck*/] = {{117, 98},  {133, 18},  {133, 130},  {140, 123},  {147, 145},  {160, 141},  {161, 159},  {169, 6},  {175, 163},  {189, 169},  {210, 67},  {214, 39},  {217, 48},  {237, 62},  {241, 196},  {245, 209},  {252, 23}}; // cyclos 3,4,5,6 all have a large prime factor

    constexpr int SizeLastCheck = sizeof(LastCheck)/sizeof(SmallRat);

    // Assume f is non-constant with integer coeffs (ideally squarefree, content=1, and f palindromic)
    VerboseLog VERBOSE("FindCycloFactor");
    const bool palindromic = IsPalindromic(f);

    const RingElem& X = indet(owner(f), UnivariateIndetIndex(f));
    const ring ZZx = NewPolyRing(RingZZ(), symbols("x")); // used only to evaluate a cyclo poly
    const RingElem& x = indet(ZZx,0);                     //


    VERBOSE(80) <<"NUM INIT Candidates: " << len(CandidateCycloIndices) << std::endl;
    long EvalPtNumer = 2;  long EvalPtDenom = 1; /// always have: EvalPtNumer > EvalPtDenom & gcd = 1; also need EvalPtNumer^3 < MAX_LONG
///???    int NumEvalPts = 0;
    bool LastIter = false;
    while (true)
    {
      BigInt Valf = EvalAt(f, EvalPtNumer, EvalPtDenom);
      while (IsZero(Valf))  // if (IsZero(Valf)) {...}   if input is square-free
      {
        f /= EvalPtDenom*X - EvalPtNumer; // exact division!    MAYBE DIVIDE BY EvalPtDenom*EvalPtNumer*(X*X+1) - X*(power(EvalPtNumer,2)+power(EvalPtDenom,2)); ???
        Valf = EvalAt(f, EvalPtNumer, EvalPtDenom);
      }
      if (!palindromic)  // if input is not palindromic...
        Valf = gcd(Valf, EvalAt(f, EvalPtDenom, EvalPtNumer));
///???      ++NumEvalPts;
      VERBOSE(80) << "Chosen EvalPt=" << EvalPtNumer << "/" << EvalPtDenom << std::endl;
      BigInt Valf_reduced = CoprimeFactor(Valf, EvalPtNumer*(EvalPtNumer*EvalPtNumer-EvalPtDenom*EvalPtDenom)); // BUG?????  to avoid overflow need EvalPt^3 < MaxLong (or MaxULong)

      vector<long> ReducedCandidateList;
      for (long k:  CandidateCycloIndices)
      {
        if (k < 3)  continue; // skip 1 & 2, if present
        if (k > 6 && IsOne(Valf_reduced))   break;
        if (EvalPtNumer == 2 && EvalPtDenom == 1 && k == 6 && Valf%3 == 0)
        { ReducedCandidateList.push_back(k); continue; } // exception in Zsygmondi's Thm
        const BigInt g = gcd(Valf_reduced,  power(EvalPtNumer,k) - power(EvalPtDenom,k));
        if (IsOne(g))  continue;
        Valf_reduced = CoprimeFactor(Valf_reduced, g);
        if (IsDivisible(Valf, EvalAt(CyclotomicPoly(k,x),EvalPtNumer,EvalPtDenom)))
          ReducedCandidateList.push_back(k);
      }
      if (LastIter)  return ReducedCandidateList;
      if (len(ReducedCandidateList) == len(CandidateCycloIndices))
      {
        // Do one last check with an eval pt chosen so that Phi(3), Phi(4) & Phi(6) have "large" prime factors
        LastIter = true;
        const long final = RandomLong(0, SizeLastCheck-1);
        EvalPtNumer = LastCheck[final].num;
        EvalPtDenom = LastCheck[final].den;
        continue;
      }
      // We filtered out at least one index...
      swap(CandidateCycloIndices, ReducedCandidateList); // really assignment
      VERBOSE(80) << "CandidateCycloIndices=" << CandidateCycloIndices << std::endl;
      if (CandidateCycloIndices.empty())  return CandidateCycloIndices;
      // advance to next EvalPt
      const long skip = RandomLong(1, 50); for (long i=0; i < skip; ++i)  GotoNextEvalPt(EvalPtNumer, EvalPtDenom);
    }
  }



  namespace // anonymous -- auxiliaries for CyclotomicFactorIndices
  {

    // List of all k>2 such that EulerTotient(k) <= d
    // EQUIV list of all indices k>2 such that deg(cyclo(k)) <= d
    std::vector<long> CycloIndicesUptoDeg(long d)
    {
      vector<long> CandidateIndices;
      const long UPB = InvTotientBoundUpto_ulong(d);
      for (long k=3; k <= UPB; ++k)
      {
        if (IsOdd(k) && k > UPB/2)  continue;
        const long phi = EulerTotient(k);
        if (phi <= d)
          CandidateIndices.push_back(k);
      }
      return CandidateIndices;
    }

  } // end of namespace anonymous


  // NOTE: this produces an UNVERIFIED list of candidate indices (but false positives are "rare")
  std::vector<long> CyclotomicFactorIndices(RingElem f)
  {
    const char* const FnName = "CyclotomicFactorIndices";
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    vector<long> IndexList; // will contain result
    if (IsConstant(f))  return IndexList; // must be before call to UnivariateIndetIndex
    const long x_index = UnivariateIndetIndex(f);
    if (x_index < 0)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    VerboseLog VERBOSE(FnName);
    VERBOSE(80) << "Minor preprocessing" << std::endl;
    // TODO: if ConstCoeff == 0 then divide by power of x
    f = prim(f);
    RingElem revf = reverse(f); // may have lower degree than f
    f = (deg(revf) < 100) ? gcd(f,revf) : revf;
    if (IsConstant(f))  return IndexList;
    const bool palindromic = (deg(f) <= deg(revf) || IsPalindromic(f));
    const RingElem& x = indet(owner(f), x_index);
    // It is convenient to deal with Phi_1 and Phi_2 separately
    if (IsZero(EvalAt(f,1)))   { f /= x-1; IndexList.push_back(1); }
    if (IsZero(EvalAt(f,-1)))  { f /= x+1; IndexList.push_back(2); }
    if (deg(f) <= 1)  return IndexList;
    VERBOSE(80) << "Checking if it is cyclotomic" << std::endl;
    const long k = CyclotomicTest(f);
    if (k != 0)  { IndexList.push_back(k);  return IndexList; }

    // Possibly reduce max poss deg of any prod of cyclo factors:
    BigInt G; // init zero
    BigInt Gnew;
    long deg_drop1 = 0; // will be set in 1st loop below
    long deg_drop2 = 0; // will be set in 2nd loop below
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      Gnew = gcd(G, ConvertTo<BigInt>(coeff(it)));
      if (IsOne(Gnew))  { deg_drop1 = deg(f) - deg(PP(it));  break; }
      swap(G, Gnew); // equiv G = Gnew;
    }
    if (palindromic)  { deg_drop2 = deg_drop1; }
    else
    {
      for (SparsePolyIter it=BeginIter(revf); !IsEnded(it); ++it)
      {
        Gnew = gcd(G, ConvertTo<BigInt>(coeff(it)));
        if (IsOne(Gnew))  { deg_drop2 = deg(revf) - deg(PP(it));  break; }
        swap(G, Gnew); // equiv G = Gnew;
      }
    }
    VERBOSE(80) << "Degree drop is " << deg_drop1+deg_drop2 << std::endl;
    if (deg(f) - deg_drop1 - deg_drop2 < 2)  return IndexList;

    vector<long> indices = FindCycloFactor(CycloIndicesUptoDeg(deg(f)-deg_drop1-deg_drop2), f); // recall that coeffs of f are integer!
    IndexList.insert(IndexList.end(), indices.begin(), indices.end()); // concat
    sort(IndexList.begin(), IndexList.end());
    return IndexList;
  }


  factorization<RingElem> CyclotomicFactors(ConstRefRingElem f)
  {
    const char* const FnName = "CyclotomicFactors";
    VerboseLog VERBOSE(FnName);
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    const factorization<RingElem> facs = SqFreeFactor(PalindromicFactor(f));
    const vector<RingElem>& F = facs.myFactors();
    const vector<long>& m = facs.myMultiplicities();
    const ring& P = owner(f);
    const RingElem& x = indet(P, UnivariateIndetIndex(f));
    factorization<RingElem> ans(one(P)); // remaining factor
    for (int i=0; i < len(F); ++i)
    {
      VERBOSE(85) << "LOOP("<<i<<"): doing " << LPP(F[i]) << " + ..." << std::endl;
      // ???MAYBE??? long k = CyclotomicIndex(F[i]); if (k != 0) { ans.myAppend(F[i], m[i]); continue; }
      const vector<long> CandidateIndices = CyclotomicFactorIndices(F[i]);
      // MAY GIVE FALSE POSITIVE { long sum = 0; for (long k: CandidateIndices) sum += EulerTotient(k); if (sum == deg(F[i])) { ...???... }}
      for (long k: CandidateIndices)
      {
        const RingElem phi_k = CyclotomicPoly(k,x);
        if (!IsDivisible(F[i], phi_k)) continue;
        ans.myAppend(phi_k, m[i]);
      }
    }
    return ans;
  }



  // -------------------------------------------------------
  // CyclotomicIndex

  namespace  // anonymous
  {

    RingElem RootX(RingElem f, long n) /*noexcept*/
    {
      // ASSUMES input is non constant, univariate.
      // ASSUMES input is of form g(x^n)
      CoCoA_ASSERT(n > 1);
      const ring& P = owner(f);
      CoCoA_ASSERT(IsPolyRing(P));
      CoCoA_ASSERT(!IsConstant(f));
      CoCoA_ASSERT(UnivariateIndetIndex(f) >= 0);

      PPMonoidElem x = radical(LPP(f));
      RingElem ret = zero(P);
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        const long d = deg(PP(it));
        ret += monomial(P, coeff(it), power(x,d/n)); // division ASSUMED exact!
      }
      return ret;
    }

  } // end of namespace anonymous



  // If f is monic cyclo with index n, return n; o/w return 0.
  // Uses CyclotomicIndex_unchecked (above) to do most of the work
  unsigned long CyclotomicIndex(ConstRefRingElem f)
  {
    return CyclotomicTest(f, false/*no full check*/);
  }

  unsigned long CyclotomicTest(ConstRefRingElem f)
  {
    return CyclotomicTest(f, true/*do full check*/);
  }


}  // end of namespace CoCoA
