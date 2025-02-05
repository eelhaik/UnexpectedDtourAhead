#define UNICODE // Enable Unicode support

#include <windows.h>
#include <shobjidl.h>   // For SHBrowseForFolder
#include <shlobj.h>      // For SHBrowseForFolder
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <cstring>  // or <string.h>
#include <cmath>
#include <iomanip>

#ifdef _MSC_VER
#define safe_strcpy(destination, size, source) strcpy_s(destination, size, source)
#else
#define safe_strcpy(destination, size, source) strncpy(destination, source, size - 1); destination[size - 1] = '\0'
#endif


// This script calculates a veriety overall measures relating to D, partitioned by transitions and transversions and by the mutating triplet.
// Output include a full table of both site counts and D values (as ABBA and BABA counts) in the upper half.
// Below the full output table are generated specific outputs that correspond to Tables / plots found in the text.
// The user is able to select which populations / regions to analyse and which (contiguous) chromosomes to analyse.

using namespace std;

std::wstring SelectFolder();

void LoadPopsHGDP(std::wstring SF);
void LoadPopDetailsHGDP(std::wstring SF);
void LoadReference(std::wstring SF, int chrom);
void Window_averages(ofstream& out);
float Correlation();
int triplet(long loc, char BS);
int base_to_num(char B);
int TV(char A, char B);
void Spike_rotate(int geno[1000][2], char ref, char alt, int trip);
void Branch_lengths(int geno[1000][2], char ref, char alt, int trip);
void Get_Run_parameters();
void Generate_output(ofstream& out);

char REF[260000000];
char popnams[100][50];
int pops[1000][2]{}; // slight array size overkill for HGDP
char line[100000];
char triptype[32][5];
int t1, t2, P1, P2, st, nd, reps = 0;
long window = 50000;
float spike = 0.02;
long types[128][11]{};
int D[2]{};
int BCH[12]{};
int Branch[3][2]{};
float ABBA[20]{};
float BLEN[3][2]{};
float Het[60]{};
float ABBAtype[500][64][12]{};
double DPrime[2000][10][2]{};
float BCH_byCount[100][7][2]{};
float DA[54]{};
int REG[54]{};
float Slope_v_D[200][6]{};


int main() {
    // Get folder path from the user
    std::string outFilename;
    std::wstring selectedFolder = SelectFolder();
    if (selectedFolder.back() != L'\\' && selectedFolder.back() != L'/') selectedFolder += L'\\';  // Add backslash if needed

    LoadPopsHGDP(selectedFolder); // load population information, including the population and region identities of all samples
    LoadPopDetailsHGDP(selectedFolder); // load details of the populations, including distance from Africa
    Get_Run_parameters(); // load run paramters

    long locHUM = 0, ct = 0; // variables for locations of each site
    char alt[10000], ref[10000]; // big character arrays to take maximum allele sequences!
    int seed = 0, prev = 0, lcn = 0, prev2 = 0, lcn2 = 0;
    char temp[10000];

    std::cout << "Random number seed "; // generate run-specific seed
    std::cin >> seed;
    for (int i = 0; i < seed; i++) lcn = rand();
std:cout << "Output file name "; // determine output filename and add '.txt' if omitted
    std::cin >> outFilename;

    // open file for output and check that it opens
    std::string outfile(selectedFolder.begin(), selectedFolder.end());
    outfile += outFilename;
    outfile += ".txt";
    std::ofstream out(outfile);
    if (out) std::cout << outfile << "  open\n";
    else std::cout << outfile << " is not open\n", exit(1); // if not open, stall

    for (int chrom = st; chrom <= nd; chrom++) {
        prev = 0, prev = 0;
        LoadReference(selectedFolder, chrom); // load full reference sequence into array REF
        std::string CR = std::to_string(chrom); // convert chromosome number to string

        // build name of Bergstrom HGDP file, open and read lines until data are reached
        std::string infile(selectedFolder.begin(), selectedFolder.end());
        infile += "hgdp.v0.5.archaics.chr" + CR + ".vcf";

        ifstream in(infile);
        if (in) cout << infile << " is open\n";
        else cout << infile << " is not open, exiting program, make sure this file is in the active folder\n", exit(2);
        while (in.peek() == '#') in.getline(line, 1000000);
        cout << "\n\n";

        while (!in.eof()) {
            in >> temp >> locHUM >> temp >> ref >> alt >> temp >> temp >> temp >> temp;
            in.getline(line, 1000000);

            lcn = locHUM / window;
            lcn2 = locHUM / 1000000;
            if (prev2 < lcn2) cout << "chomosome " << chrom << ", megabases " << prev2 << "\n";
            if (prev < lcn) Window_averages(out);

            prev = lcn;
            prev2 = lcn2;

            if (temp[0] == 'G' && temp[1] == 'T') { // check that 'line' contains gnotype data
                int na = 0, at = 0; // number of alleles scored and number of alternate alleles
                int geno[1000][2]{}; // genotypes stored in array
                for (int j = 0; j < 934; j++) {
                    if (line[j * 4 + 1] != '.' && line[j * 4 + 3] != '.') { // only translate scored individuals
                        geno[j][0] = line[j * 4 + 1] - 48; // convert to '0' and '1'
                        geno[j][1] = line[j * 4 + 3] - 48;
                        if (j > 4) at += (geno[j][0] + geno[j][1]); // count alternate alleles
                        na += 2; // count scored alleles
                    }
                    else geno[j][0] = -1, geno[j][1] = -1; // if not scored, allele stored as -1
                }

                // check first 5 genotypes are all homozygous and scored
                bool test = true;
                for (int i = 0; i < 5; i++) if ((geno[i][0] != geno[i][1]) || (geno[i][0] == -1)) test = false;

                if (na > 1800 && test && ref[1] == '\0' && alt[1] == '\0') { // check biallelic, more than 900 scored, non-humans homozygous, some alternates
                    int TP = 32;
                    if (at < float(na) / 2) TP = triplet(locHUM, ref[0]); // determine triplet of current site based on major allele mutating to minor
                    else TP = triplet(locHUM, alt[0]);
                    int ts = TV(ref[0], alt[0]); // determine transition to transversion
                    TP = TP + 32 * ts;
                    Spike_rotate(geno, ref[0], alt[0], TP); // analyse the chosen population combination for D classified by type
                    Branch_lengths(geno, ref[0], alt[0], TP);
                    ct++;
                    if (ct >= 100000) reps++, ct = 0; // generate blocks of 100000 contiguous qualifying sites
                }
            }
        }
    }

    Generate_output(out);
    return 0;
}

void Window_averages(ofstream& out)
{
    if (BCH[0] > 2500) {
        int cdist = (25 * (float(BCH[0] - BCH[1]) / (BCH[0] + BCH[1]) + 1));
        DPrime[cdist][0][0] += BCH[2], DPrime[cdist][0][1] += BCH[3]; // dprime 0
        cdist = (500 * (float(BCH[0] - BCH[1]) / (BCH[0] + BCH[1]) + 1));
        DPrime[cdist][6][0] += BCH[2], DPrime[cdist][6][1] += BCH[3]; // dprime 6

        cdist = (25 * (float(BCH[4] - BCH[5]) / (BCH[4] + BCH[5]) + 1));
        DPrime[cdist][1][0] += BCH[6], DPrime[cdist][1][1] += BCH[7]; // dprime 1
        cdist = (500 * (float(BCH[4] - BCH[5]) / (BCH[4] + BCH[5]) + 1));
        DPrime[cdist][7][0] += BCH[6], DPrime[cdist][7][1] += BCH[7]; // dprime 7

        cdist = (25 * (float(BCH[8] - BCH[9]) / (BCH[8] + BCH[9]) + 1));
        DPrime[cdist][2][0] += BCH[10], DPrime[cdist][2][1] += BCH[11]; // dprime 2

        DPrime[cdist][3][0] += ABBA[0], DPrime[cdist][3][1] += ABBA[1];
        DPrime[cdist][4][0] += Het[58], DPrime[cdist][4][1] += Het[59], Het[58] = 0, Het[59] = 0;
        DPrime[cdist][5][0] += D[0], DPrime[cdist][5][1] += D[1];

        cdist = (500 * (float(BCH[8] - BCH[9]) / (BCH[8] + BCH[9]) + 1));
        DPrime[cdist][8][0] += BCH[10], DPrime[cdist][8][1] += BCH[11]; // dprime 8
    }
    if (Branch[0][0] + Branch[0][1] > 0) {
        int DCP1 = Branch[0][0];
        int DCP2 = Branch[0][1];
        int cdist = (25 * (float(DCP1 - DCP2) / (DCP1 + DCP2) + 1));
        for (int i = 0; i < 3; i++) {
            BCH_byCount[cdist][i][0] += Branch[i][0];
            BCH_byCount[cdist][i][1] += Branch[i][1];
        }
        BCH_byCount[cdist][3][0] += ABBA[0];
        BCH_byCount[cdist][3][1] += ABBA[1];
    }
    if (Het[6] + Het[43] > 0) {
        int asym = (50 * (float(Het[6] - Het[43]) / (Het[6] + Het[43]) + 1));
        for (int i = 0; i < 3; i++) {
            BCH_byCount[asym][i+4][0] += BLEN[i][0];
            BCH_byCount[asym][i+4][1] += BLEN[i][1];
        }
    }

    int slope = (10000 * Correlation()) + 50;
    if (slope < 0) slope = 0;
    if (slope > 99) slope = 99;
    int Hetdif = Het[0] - Het[27] + 10;
    if (Hetdif < 0) Hetdif = 0;
    if (Hetdif > 89) Hetdif = 89;
    for (int g = 0; g < 6; g++) Slope_v_D[Hetdif + 100][g] += ABBA[g];
    for (int g = 0; g < 6; g++) Slope_v_D[slope][g] += ABBA[g], ABBA[g] = 0;

    for (int i = 0; i < 12; i++) BCH[i] = 0;
    D[0] = 0, D[1] = 0;
    for (int i = 0; i < 3; i++) Branch[i][0] = 0, Branch[i][1] = 0, BLEN[i][0] = 0, BLEN[i][1] = 0;
    for (int i = 0; i < 54; i++) Het[i] = 0;
    // out << "\n";
}

int TV(char A, char B)
{
    if ((A == 'A' && B == 'G') || (A == 'G' && B == 'A')) return 1;
    if ((A == 'C' && B == 'T') || (A == 'T' && B == 'C')) return 1;
    return 0;
}

float Correlation()
{
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
    int n = 0;

    // Calculate the necessary sums
    for (int i = 0; i < 54; ++i) {
        if (REG[i] == 0 || REG[i] == 1 || REG[i] == 3 || REG[i] == 4) {
            sumX += DA[i];
            sumY += Het[i];
            sumXY += DA[i] * Het[i];
            sumX2 += DA[i] * DA[i];
            n++;
        }
    }

    // Calculate the slope (m) using the formula for linear regression
    double numerator = n * sumXY - sumX * sumY;
    double denominator = n * sumX2 - sumX * sumX;

    if (denominator == 0) {
        std::cerr << "Error: Denominator is zero, unable to compute slope." << std::endl;
        return 0.0;
    }

    double slope = numerator / denominator;
    return slope;
}

int base_to_num(char B)
{
    if (B == 'A') return 0;
    else if (B == 'C') return 1;
    else if (B == 'G') return 2;
    else if (B == 'T') return 3;
    else return -1;
}

int triplet(long loc, char BS) // determine triplet as a number, converting triplets with central G or T to opposite strand A or C
{
    char b1 = base_to_num(REF[loc - 1]);
    char b2 = base_to_num(REF[loc + 1]);
    if (b1 >= 0 && b2 >= 0) {
        if (BS == 'A') return 4 * b1 + b2;
        else if (BS == 'C') return 16 + 4 * b1 + b2;
        else if (BS == 'G') return 16 + 4 * (3 - b2) + (3 - b1);
        else if (BS == 'T') return 4 * (3 - b2) + (3 - b1);
        else return 32;
    }
    else return 32;
}

void Branch_lengths(int geno[1000][2], char ref, char alt, int trip) // calculate branch lengths, heterozygositys and ABBAs / BABAs for a given informative site
{
    float freq[54][2]{};
    float FQ[7][2]{};

    int a = 0, n = 0;
    for (int i = 5; i < 934; i++) { // determine population (freq[pop][]) and regional (FQ[region][]) alternate allele counts [0] = aa, [1] = all alleles
        if (geno[i][0] == 1) freq[pops[i][0]][0]++, FQ[pops[i][1]][0]++, a++;
        if (geno[i][1] == 1) freq[pops[i][0]][0]++, FQ[pops[i][1]][0]++, a++;
        if (geno[i][0] != -1) freq[pops[i][0]][1] += 2, FQ[pops[i][1]][1] += 2, n += 2; // only count scored alleles
    }

    if ((t1 == 0 && freq[P1][1] == 0) || (t2 == 0 && freq[P2][1] == 0)) return;

    float fP1, fP2;
    if (t1 == 0) fP1 = float(freq[P1][0]) / freq[P1][1];
    else fP1 = float(FQ[P1][0]) / FQ[P1][1];
    if (t2 == 0) fP2 = float(freq[P2][0]) / freq[P2][1];
    else fP2 = float(FQ[P2][0]) / FQ[P2][1];

    if (fP1 == 0 && fP2 > 0) { // AFR = ref, nAFR carries some alt
        if (geno[0][0] == 0 && geno[3][0] == 0) Branch[0][1]++, BLEN[0][1] += fP2; // P2 has mutated relative to cmp & nea
        if (geno[0][0] == 0 && geno[3][0] == 1) Branch[1][1]++, BLEN[1][1] += fP2; // P2 has mutated relative to cmp but not nea
        if (geno[0][0] == 1 && geno[3][0] == 0) Branch[2][1]++, BLEN[2][1] += fP2; // P2 has mutated relative to nea but not cmp
    }
    if (fP1 > 0 && fP2 == 0) {
        if (geno[0][0] == 0 && geno[3][0] == 0) Branch[0][0]++, BLEN[0][0] += fP1; // P1 has mutated relative to cmp & nea
        if (geno[0][0] == 0 && geno[3][0] == 1) Branch[1][0]++, BLEN[1][0] += fP1; // P1 has mutated relative to cmp but not nea
        if (geno[0][0] == 1 && geno[3][0] == 0) Branch[2][0]++, BLEN[2][0] += fP1; // P1 has mutated relative to nea but not cmp
    }
    if (fP1 == 1 && fP2 < 1) {
        if (geno[0][0] == 1 && geno[3][0] == 1) Branch[0][1]++, BLEN[0][1] += (1 - fP2);
        if (geno[0][0] == 1 && geno[3][0] == 0) Branch[1][1]++, BLEN[1][1] += (1 - fP2);
        if (geno[0][0] == 0 && geno[3][0] == 1) Branch[2][1]++, BLEN[2][1] += (1 - fP2);
    }
    if (fP1 < 1 && fP2 == 1) {
        if (geno[0][0] == 1 && geno[3][0] == 1) Branch[0][0]++, BLEN[0][0] += (1 - fP1);
        if (geno[0][0] == 1 && geno[3][0] == 0) Branch[1][0]++, BLEN[1][0] += (1 - fP1);
        if (geno[0][0] == 0 && geno[3][0] == 1) Branch[2][0]++, BLEN[2][0] += (1 - fP1);
    }

    if (geno[0][0] == 0 && geno[3][0] == 1) {
        ABBA[0] += (1 - fP1) * fP2; // P2 = Neanderthal
        ABBA[1] += fP1 * (1 - fP2); // P1 = Neanderthal
        if (a < float(n) / 2) ABBA[2] += (1 - fP1) * fP2, ABBA[3] += fP1 * (1 - fP2); // major allele = chimpanzee allele
        else ABBA[4] += (1 - fP1) * fP2, ABBA[5] += fP1 * (1 - fP2); // major allele = Neanderthal allele
    }
    if (geno[0][0] == 1 && geno[3][0] == 0) {
        ABBA[0] += fP1 * (1 - fP2); // P2 = Neanderthal
        ABBA[1] += (1 - fP1) * fP2; // P1 = Neanderthal
        if (a > float(n) / 2) ABBA[2] += fP1 * (1 - fP2), ABBA[3] += (1 - fP1) * fP2; // major allele = chimpanzee allele
        else ABBA[4] += fP1 * (1 - fP2), ABBA[5] += (1 - fP1) * fP2; // major allele = Neanderthal allele
    }

    Het[58] += 2 * fP1 * (1 - fP1);
    Het[59] += 2 * fP2 * (1 - fP2);

    for (int pop = 0; pop < 54; pop++) {
        if (freq[pop][1] > 0) {
            float fq = float(freq[pop][0]) / freq[pop][1];
            Het[pop] += 2 * fq * (1 - fq);
        }
    }
}

void Spike_rotate(int geno[1000][2], char rf, char at, int trip) // generate ABBA and BABA counts and branch lengths for each 'block' (defined by 'reps')
{
    int cnt1 = 0, arch = 0;
    if (geno[3][0] == geno[2][0] && geno[3][0] == geno[4][0]) arch = 6;
    for (int i = 0; i < 929; i++) {
        if ((t1 == 0 && pops[i + 5][0] == P1) || (t1 == 1 && pops[i + 5][1] == P1)) {
            int cnt2 = 0, cnt3 = 0;
            for (int j = 0; j < 929; j++) {
                bool test1 = ((t1 == 0 && pops[j + 5][0] == P1) || (t1 == 1 && pops[j + 5][1] == P1));
                bool test2 = ((t2 == 0 && pops[j + 5][0] == P2) || (t2 == 1 && pops[j + 5][1] == P2));
                if ((test1 || test2) && j != i) {
                    int A1, A2, A3;
                    if (float(rand()) / RAND_MAX < 0.5) A1 = geno[i + 5][0];
                    else A1 = geno[i + 5][1];
                    if (float(rand()) / RAND_MAX < 0.5) A2 = geno[j + 5][0];
                    else A2 = geno[j + 5][1];
                    if (float(rand()) / RAND_MAX < spike) A3 = geno[3][0];
                    else A3 = A2;

                    if (pops[i + 5][0] == pops[j + 5][0]) {  // sampe population, conduct spiked and unspiked analysis                       
                        if (geno[0][0] == A1 && geno[0][0] != A2) BCH[1]++; //branch to cmp P1 v P1
                        if (geno[0][0] == A2 && geno[0][0] != A1) BCH[0]++;
                        if (geno[0][0] == A1 && geno[0][0] != A3) BCH[5]++; //branch to cmp P1 v P1 spiked
                        if (geno[0][0] == A3 && geno[0][0] != A1) BCH[4]++;

                        if (geno[3][0] == A1 && geno[3][0] != A2) BCH[3]++; //branch to Nea P1 v P1 
                        if (geno[3][0] == A2 && geno[3][0] != A1) BCH[2]++;
                        if (geno[3][0] == A1 && geno[3][0] != A3) BCH[7]++; //branch to Nea P1 v P1 spiked
                        if (geno[3][0] == A3 && geno[3][0] != A1) BCH[6]++;

                        if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 0 && A2 == 1) ABBAtype[reps][trip][2 + arch]++;
                        if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 1 && A2 == 0) ABBAtype[reps][trip][3 + arch]++;
                        if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 0 && A2 == 1) ABBAtype[reps][trip][3 + arch]++;
                        if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 1 && A2 == 0) ABBAtype[reps][trip][2 + arch]++;

                        if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 0 && A3 == 1) ABBAtype[reps][trip][4 + arch]++;
                        if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 1 && A3 == 0) ABBAtype[reps][trip][5 + arch]++;
                        if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 0 && A3 == 1) ABBAtype[reps][trip][5 + arch]++;
                        if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 1 && A3 == 0) ABBAtype[reps][trip][4 + arch]++;
                    }
                    else {
                        if (test2) {
                            if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 0 && A2 == 1) D[1]++, ABBAtype[reps][trip][0 + arch]++;
                            if (geno[0][0] == 0 && geno[3][0] == 1 && A1 == 1 && A2 == 0) D[0]++, ABBAtype[reps][trip][1 + arch]++;
                            if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 0 && A2 == 1) D[0]++, ABBAtype[reps][trip][1 + arch]++;
                            if (geno[0][0] == 1 && geno[3][0] == 0 && A1 == 1 && A2 == 0) D[1]++, ABBAtype[reps][trip][0 + arch]++;

                            if (geno[0][0] == A1 && geno[0][0] != A2) BCH[9]++; //branch to cmp from P2
                            if (geno[0][0] == A2 && geno[0][0] != A1) BCH[8]++; //branch to cmp from P1
                            if (geno[3][0] == A1 && geno[3][0] != A2) BCH[11]++; //branch to Nea P1 v P2
                            if (geno[3][0] == A2 && geno[3][0] != A1) BCH[10]++;
                        }
                    }
                }
            }
            cnt1++;
        }
    }
}

void LoadReference(std::wstring SF, int chrom) // loads reference sequence and counts the numbers of each triplet (TRIPS[])
{
    std::string CR = std::to_string(chrom);
    std::string infile(SF.begin(), SF.end());
    infile += "chr" + CR + ".fa";
    std::ifstream in2(infile);
    cout << infile << "\n";
    if (!in2) cout << "not open " << infile << "\n";
    cout << "Loading the reference sequence . . \n";
    for (long i = 0; i < 260000000; i++) REF[i] = 'N'; // all bases are 'N' unless scored

    int cr = 0;
    in2.getline(line, 1000);

    long cnt = 0;
    long c = 0;
    char base;
    // while (!in2.eof() && cnt < 249999999) {
    while (!in2.eof()) {
        base = in2.get();
        if (int(base) > 64 && int(base) < 87) {
            cnt++;
            c++;
            if (c > 9999999) cout << cnt << "\n", c = 0;
            REF[cnt] = base;
        }
        else if (int(base) > 96 && int(base) < 120) {
            char typ = 'N';
            if (base == 'a') typ = 'A';
            else if (base == 'c') typ = 'C';
            else if (base == 'g') typ = 'G';
            else if (base == 't') typ = 'T';
            cnt++;
            c++;
            if (c > 9999999) cout << cnt << "\n", c = 0;
            REF[cnt] = typ;
        }
    }
    in2.close();
    in2.clear();
}

void LoadPopDetailsHGDP(std::wstring SF)
{
    std::string infile(SF.begin(), SF.end());
    infile += "HGDP_pop_characters.txt";
    std::ifstream in1(infile);
    if (!in1) cout << infile << "  NOT OPEN\n";
    in1.getline(line, 1000);
    char temp[100];

    for (int i = 0; i < 54; i++) {
        in1 >> temp >> temp >> temp >> DA[i] >> REG[i];
    }
    in1.close();
    in1.clear();
}

void LoadPopsHGDP(std::wstring SF) // read in population and region codes for each of 929 individuals
{
    int N_ind[60]{};
    int N_reg[60]{};
    for (int i = 0; i < 100; i++) popnams[i][0] = 'w';
    std::string infile(SF.begin(), SF.end());
    infile += "inpopsHGDP.txt";
    std::ifstream in1(infile);

    char temp[100], nam[100];
    for (int i = 0; i < 929; i++) {
        in1 >> temp >> nam >> temp >> temp >> pops[i + 5][0] >> pops[i + 5][1]; // first 5 individuals are not HGDP samples, hence +5
        pops[i + 5][0]--; // population codes start at 1 so decrement for use as array indices
        pops[i + 5][1]--;
        N_ind[pops[i + 5][0]]++;
        N_reg[pops[i + 5][1]]++;
        if (popnams[pops[i + 5][0]][0] == 'w') safe_strcpy(popnams[pops[i + 5][0]], sizeof(popnams[pops[i + 5][0]]), nam);
    }
    cout << "popnum" << "\t" << "name" << "\t\t" << "N" << "\n";
    for (int pop = 0; pop < 54; pop += 3) {
        for (int i = 0; i < 3; i++) {
            std::cout << pop + i << "\t" << popnams[pop + i] << "\t";
            if (strlen(popnams[pop + i]) < 8) cout << "\t";
            std::cout << N_ind[pop + i];
            if (i < 2) cout << "\t";
        }
        std::cout << "\n";
    }

    // store region and population names
    cout << "\n";
    std::cout << 0 << "\t" << "Africa" << "\t" << N_reg[0] << "\t" << 1 << "\t" << "MidEast" << "\t" << N_reg[1] << "\t" << 2 << "\t" << "Europe" << "\t" << N_reg[2] << "\n";
    std::cout << 3 << "\t" << "S Asia" << "\t" << N_reg[3] << "\t" << 4 << "\t" << "E Asia" << "\t" << N_reg[4] << "\t" << 5 << "\t" << "Oceania" << "\t" << N_reg[5] << "\n";
    std::cout << 6 << "\t" << "America" << "\t" << N_reg[6] << "\n";

    safe_strcpy(popnams[60], sizeof(popnams[60]), "Africa");
    safe_strcpy(popnams[61], sizeof(popnams[61]), "MidEast");
    safe_strcpy(popnams[62], sizeof(popnams[62]), "Europe");
    safe_strcpy(popnams[63], sizeof(popnams[63]), "S Asia");
    safe_strcpy(popnams[64], sizeof(popnams[64]), "E Asia");
    safe_strcpy(popnams[65], sizeof(popnams[65]), "Oceania");
    safe_strcpy(popnams[66], sizeof(popnams[66]), "America");
    safe_strcpy(popnams[70], sizeof(popnams[70]), "Vindija");
    safe_strcpy(popnams[71], sizeof(popnams[71]), "Altai");
    safe_strcpy(popnams[72], sizeof(popnams[72]), "Denisova");
    in1.close();
}

void Get_Run_parameters()
{
    t1 = 0, P1 = 6, t2 = 0, P2 = 43, st = 22, nd = 22;
    char k = 'z';
    cout << "\nInput run parameters by choosing the populations / regions to analyse \n";
    while (k != 'Y') {
        cout << "P1 population or region\t\t" << t1 << "\t";
        if (t1 == 0) cout << "popn" << "\t" << "\t<a>\n";
        else cout << "region" << "\t" << "\t<a>\n";
        cout << "P1 taxon code\t\t\t" << P1 << "\t" << popnams[P1 + t1 * 60] << "\t" << "\t<b>\n";
        cout << "P2 population or region\t\t" << t2 << "\t";
        if (t2 == 0) cout << "popn" << "\t" << "\t<c>\n";
        else cout << "region" << "\t" << "\t<c>\n";
        cout << "P2 taxon code\t\t\t" << P2 << "\t" << popnams[P2 + t2 * 60] << "\t" << "\t<d>\n";
        cout << "start chromosome\t\t" << st << "\t\t" << "\t<e>\n";
        cout << "end chromosome\t\t\t" << nd << "\t\t" << "\t<f>\n";
        cout << "window size\t\t\t" << window << "\t\t\t<g>\n";
        cout << "spiking percent\t\t\t" << spike << "\t\t\t<h>\n\n";
        cout << "\t\tto end entry and run <Y> to terminate program <Z>  ";
        cin >> k;
        if (k == 'Z') exit(0);
        if (k == 'a') t1 = 1 - t1;
        if (k == 'c') t2 = 1 - t2;
        if (k == 'a' || k == 'b') {
            cout << "enter new taxon code for P1 ";
            cin >> P1;
        }
        if (k == 'c' || k == 'd') {
            cout << "enter new taxon code for P2 ";
            cin >> P2;
        }
        if (t1 == 1 && P1 > 6) {
            cout << "\nP1 is a region so code must be between 0 and 6, please re-enter\n";
            cin >> P1;
        }
        if (t2 == 1 && P2 > 6) {
            cout << "\nP2 is a region so code must be between 0 and 6, please re-enter\n";
            cin >> P2;
        }
        if (k == 'e') {
            cout << "\nenter new start chromosome  ";
            cin >> st;
        }
        if (k == 'f') {
            cout << "\nenter new end chromosome  ";
            cin >> nd;
        }
        if (k == 'g') {
            cout << "\nenter window size  ";
            cin >> window;
        }
        if (k == 'h') {
            cout << "\nenter spiking proportion  ";
            cin >> spike;
        }

        if (k != 'a' && k != 'b' && k != 'c' && k != 'd' && k != 'e' && k != 'f' && k != 'g' && k != 'h' && k != 'Y') {
            cout << "invalid code, please re-enter ";
            cin >> k;
        }
        cout << "\n";
    }
    for (int f2 = 0; f2 < 16; f2++) triptype[f2][1] = 'A', triptype[f2 + 16][1] = 'C';
    for (int f1 = 0; f1 < 4; f1++) {
        char b;
        if (f1 == 0) b = 'A';
        if (f1 == 1) b = 'C';
        if (f1 == 2) b = 'G';
        if (f1 == 3) b = 'T';
        for (int f3 = 0; f3 < 8; f3++) triptype[f3 * 4 + f1][2] = b;
        for (int f3 = 0; f3 < 4; f3++) triptype[f1 * 4 + f3][0] = b, triptype[f1 * 4 + 16 + f3][0] = b;
    }
}

void Generate_output(ofstream& out)
{
    // calculate means and sums of squares
    float means[64][12][3]{};
    for (int trip = 0; trip < 64; trip++) {
        float D = 0;
        for (int i = 0; i < reps; i++) {
            for (int typ = 0; typ < 6; typ++) {
                if (ABBAtype[i][trip][typ * 2] + ABBAtype[i][trip][typ * 2 + 1] > 0) {
                    D = float(ABBAtype[i][trip][typ * 2] - ABBAtype[i][trip][typ * 2 + 1]) / (ABBAtype[i][trip][typ * 2 + 1] + ABBAtype[i][trip][typ * 2]);
                    means[trip][typ][0] += D;
                    means[trip][typ][1] += D * D;
                    means[trip][typ][2] ++;
                }
            }
        }
    }

    out << "Branch length and D analysis in individual windows\n\n";
    out << "P1 = ";
    if (t1 == 0) out << popnams[P1] << "\n";
    else out << popnams[P2] << "\n";
    out << "P2 = ";
    if (t2 == 0) out << popnams[P2] << "\n";
    else out << popnams[P2] << "\n";
    out << "Hetdif = Bantu - Japanese\n";
    out << "spike = " << spike << "\n";
    out << "chromosomes " << st << " to " << nd << "\n";
    out << "window size = " << window << "\n\n";

    out << "Figure 4. Branch length to Neanderthals related to branch length to chimpanzees based on locus counts\n\n";

    out << "C_dist\t" << "N>P1\t" << "N>P2\t" << "D'N\t" << "ABBA\t" << "BABA\t" << "D\n";
    for (int i = 0; i < 100; i++) {
        out << std::fixed << std::setprecision(2) << float(i)/25 - 1;
        long b1 = BCH_byCount[i][2][0], b2 = BCH_byCount[i][2][1];
        out << "\t" << b1 << "\t" << b2;
        if (b1 + b2 > 0) out << "\t" << float(b1 - b2) / (b1 + b2);
        else out << "\t";
        out << "\t" << BCH_byCount[i][3][0] << "\t" << BCH_byCount[i][3][1];
        if (BCH_byCount[i][3][0] + BCH_byCount[i][3][1] > 0) out << "\t" << float(BCH_byCount[i][3][0] - BCH_byCount[i][3][1]) / (BCH_byCount[i][3][0] + BCH_byCount[i][3][1]);
        else out << "\t";
        out << "\n";
    }
    out << "\n\n";

    out << "Figure 5, main plot. Branch length to Neanderthals related to branch length to chimpanzees based on frequencies\n\n";

    out << "C_dist\t" << "N>P1\t" << "N>P1\t" << "N>P1\t" << "N>P1_S\t" << "N>P1\t" << "N>P2\t" << "D(P1,P1)'N\t" << "D(P1,P1S)'N" << "D(P1,P2)'N\n";
    for (int i = 0; i < 100; i++) {
        out << std::fixed << std::setprecision(3) << float(i) / 25 - 1;
        for (int h = 0; h < 3; h++) out << "\t" << DPrime[i][h][0] << "\t" << DPrime[i][h][1];
        for (int h = 0; h < 3; h++) {
            if (DPrime[i][h][0] + DPrime[i][h][1] > 0) out << "\t" << float(DPrime[i][h][0] - DPrime[i][h][1]) / (DPrime[i][h][0] + DPrime[i][h][1]);
            else out << "\t";
        }
        out << "\n";
    }
    out << "\n\n";

    out << "Figure 5, inset. Branch length to Neanderthals related to branch length to chimpanzees based on frequencies\n\n";

    out << "C_dist\t" << "N>P1\t" << "N>P1\t" << "N>P1\t" << "N>P1_S\t" << "N>P1\t" << "N>P2\t" << "D(P1,P1)'N\t" << "D(P1,P1S)'N" << "D(P1,P2)'N\n";
    for (int i = 0; i < 100; i++) {
        out << std::fixed << std::setprecision(2) << float(i + 450) / 500 - 1;
        for (int h = 6; h < 9; h++) out << "\t" << DPrime[i+450][h][0] << "\t" << DPrime[i+450][h][1];
        for (int h = 6; h < 9; h++) {
            if (DPrime[i+450][h][0] + DPrime[i+450][h][1] > 0) out << "\t" << float(DPrime[i+450][h][0] - DPrime[i+450][h][1]) / (DPrime[i+450][h][0] + DPrime[i+450][h][1]);
            else out << "\t";
        }
        out << "\n";
    }
    out << "\n\n";

    out << "Figure 6, D calculated separately based on different classes of variants\n\n";
    out << "\tArchaics polymorphic\t\t\t\t\t\t" << "Archaics monomorphic\n";
    out << "Triplet\t" << "D(P1, P2)\t" << "s.e.\t" << "D(P1, P1)\t" << "s.e.\t" << "D(P1, P1_S)\t" << "s.e.\t";
    out << "D(P1, P2)\t" << "s.e.\t" << "D(P1, P1)\t" << "s.e.\t" << "D(P1, P1_S)\t" << "s.e.\n";
    for (int trip = 0; trip < 64; trip++) {
        if (trip < 32) out << triptype[trip] << "_TV";
        else out << triptype[trip - 32] << "_TS";
        for (int i = 0; i < 6; i++) {
            float se = (means[trip][i][1] - means[trip][i][0] * means[trip][i][0] / means[trip][i][2]) / (means[trip][i][2] - 1);
            se = sqrt(se / reps);
            out << "\t" << means[trip][i][0] / means[trip][i][2] << "\t" << se;
        }
        out << "\n";
    }
    out << "\n\n";

    out << "Figure 11, how D varies with the slope of the relationship between heterozygosity and distance from Africa\n\n";
    out << "\tALL\t\tMajor=A\t\tMajor=B\n";
    out << "slope\t" << "ABBA\t" << "BABA\t" << "ABBA\t" << "BABA\t" << "ABBA\t" << "BABA\t" << "D_ALL\t" << "D_MajorA\t" << "D_MajorB\n";
    for (int i = 0; i < 100; i++) {
        out << i - 50;
        for (int j = 0; j < 6; j++) out << "\t" << Slope_v_D[i][j];
        if (Slope_v_D[i][0] + Slope_v_D[i][1] > 0) out << "\t" << (Slope_v_D[i][0] - Slope_v_D[i][1]) / (Slope_v_D[i][0] + Slope_v_D[i][1]);
        else out << "\t";
        if (Slope_v_D[i][2] + Slope_v_D[i][3] > 0) out << "\t" << (Slope_v_D[i][2] - Slope_v_D[i][3]) / (Slope_v_D[i][2] + Slope_v_D[i][3]);
        else out << "\t";
        if (Slope_v_D[i][4] + Slope_v_D[i][5] > 0) out << "\t" << (Slope_v_D[i][4] - Slope_v_D[i][5]) / (Slope_v_D[i][4] + Slope_v_D[i][5]);
        else out << "\t";
        out << "\n";
    }
    out << "\n\n";

    out << "Figure 11, how D varies with the heterozygosity difference Bantu minus Japanese\n\n";
    out << "\tALL\t\tMajor=A\t\tMajor=B\n";
    out << "Het_dif\t" << "ABBA\t" << "BABA\t" << "ABBA\t" << "BABA\t" << "ABBA\t" << "BABA\t" << "D_ALL\t" << "D_MajorA\t" << "D_MajorB\n";
    for (int i = 100; i < 200; i++) {
        out << i-100;
        for (int j = 0; j < 6; j++) out << "\t" << Slope_v_D[i][j];
        if (Slope_v_D[i][0] + Slope_v_D[i][1] > 0) out << "\t" << (Slope_v_D[i][0] - Slope_v_D[i][1]) / (Slope_v_D[i][0] + Slope_v_D[i][1]);
        else out << "\t";
        if (Slope_v_D[i][2] + Slope_v_D[i][3] > 0) out << "\t" << (Slope_v_D[i][2] - Slope_v_D[i][3]) / (Slope_v_D[i][2] + Slope_v_D[i][3]);
        else out << "\t";
        if (Slope_v_D[i][4] + Slope_v_D[i][5] > 0) out << "\t" << (Slope_v_D[i][4] - Slope_v_D[i][5]) / (Slope_v_D[i][4] + Slope_v_D[i][5]);
        else out << "\t";
        out << "\n";
    }
    out << "\n\n";

    out << "Figure S3, how branch length asymmetry varies with the heterozygosity difference Yoruba minus French\n\n";
    out << "\tC=A,N=A\t\t" << "C=A,N=B\t\t" << "C=B,N=A\n";
    out << "C_dist\t" << "N>P1\t" << "N>P2\t" << "N>P1\t" << "N>P2\t" << "N>P1\t" << "N>P2\t" << "D'N(C=A,N=A)\t" << "D'N(C=A,N=B)\t" << "D'N(C=B,N=A)\n";
    for (int i = 0; i < 100; i++) {
        out << std::fixed << std::setprecision(3) << float(i) / 50 - 1;
        for (int h = 4; h < 7; h++) out << "\t" << BCH_byCount[i][h][0] << "\t" << BCH_byCount[i][h][1];
        for (int h = 4; h < 7; h++) {
            if (BCH_byCount[i][h][0] + BCH_byCount[i][h][1] > 0) out << "\t" << float(BCH_byCount[i][h][0] - BCH_byCount[i][h][1]) / (BCH_byCount[i][h][0] + BCH_byCount[i][h][1]);
            else out << "\t";
        }
        out << "\n";
    }


}


std::wstring SelectFolder() {
    // Initialize COM (required for some shell functions)
    CoInitialize(NULL);

    // Browse folder dialog
    BROWSEINFO bi;
    ZeroMemory(&bi, sizeof(bi));
    bi.lpszTitle = L"Select folder containing input files "; // Dialog title
    bi.ulFlags = BIF_NEWDIALOGSTYLE | BIF_RETURNONLYFSDIRS; // Show folders only

    // Get the selected folder path
    LPITEMIDLIST pidl = SHBrowseForFolder(&bi);
    std::wstring folderPath;

    if (pidl != NULL) {
        // Allocate memory for the path and retrieve the folder path
        wchar_t szPath[MAX_PATH];
        if (SHGetPathFromIDList(pidl, szPath)) {
            folderPath = szPath;
        }
        // Free the item ID list
        CoTaskMemFree(pidl);
    }
    else {
        std::wcout << L"No folder selected." << std::endl;
    }

    // Uninitialize COM
    CoUninitialize();

    return folderPath;
}