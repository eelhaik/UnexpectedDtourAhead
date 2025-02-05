#define UNICODE // Enable Unicode support

#include <windows.h>
#include <shobjidl.h>   // For SHBrowseForFolder
#include <shlobj.h>      // For SHBrowseForFolder
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cstring>  // or <string.h>

#ifdef _MSC_VER
#define safe_strcpy(destination, size, source) strcpy_s(destination, size, source)
#else
#define safe_strcpy(destination, size, source) strncpy(destination, source, size - 1); destination[size - 1] = '\0'
#endif

// This script calculates a variety of branch length measures for non-overlapping 50kb genomic windows, comparing African (all together) and European (all together) distances.
// The output comprises ABBA and BABA counts for each window, African - non-African branch lengths to non humans (Altai = chimpanzee), to the chimpanzee (Altai != chimpanzee) and to the Altai
// The second half of the output is based on individual Yoruba - French comparisons



using namespace std;

std::wstring SelectFolder();

void LoadPopsHGDP(std::wstring SF);
void LoadReference(std::wstring SF, int chrom);
void Get_Run_parameters();
void D_calculations(int geno[1000][2], char ref, char alt, int trip, int tstv, int major);
void Generate_output(ofstream& out);

int triplet(long loc, char BS);
int base_to_num(char B);
int TV(char A, char B);
int Get_type(int geno[1000][2], char ref, char alt, int alternate);

char REF[260000000];
char popnams[100][50];
int pops[1000][2]{}; // slight overkill for HGDP
char line[100000];
char triptype[33][5];
long types[128][12]{};
int t1, t2, P1, P2, st, nd;
float ABBA[64][36][2]{};


int main() {
    // Get folder path from the user
    std::string outFilename;
    std::wstring selectedFolder = SelectFolder();
    if (selectedFolder.back() != L'\\' && selectedFolder.back() != L'/') selectedFolder += L'\\';  // Add backslash if needed

    LoadPopsHGDP(selectedFolder);
    Get_Run_parameters();

    long locHUM = 0; // variables for locations of each site
    char alt[10000], ref[10000]; // big character arrays for max allele sequences!
    int prev = 0, prev2 = 0, lcn = 0, seed = 0, lcn2=0;
    char temp[10000];

    std:cout << "Output file name ";
    std::cin >> outFilename;

    // open file for output and check that it opens
    std::string outfile(selectedFolder.begin(), selectedFolder.end());
    outfile += outFilename;
    if (outfile.size() < 4 || outfile.substr(outfile.size() - 4) != ".txt") outfile += ".txt";
    std::ofstream out(outfile);
    if (out) std::cout << outfile << "  open\n";
    else std::cout << outfile << " is not open\n", exit(1); // if not open, stall

    for (int chrom = st; chrom <= nd; chrom++) {
        prev = 0, prev2 = 0;
        LoadReference(selectedFolder, chrom);
        std::string CR = std::to_string(chrom); // convert chromosome number to string

        // build name of Bergstrom HGDP file
        std::string infile(selectedFolder.begin(), selectedFolder.end());
        infile += "hgdp.v0.5.archaics.chr" + CR + ".vcf";

        ifstream in(infile);
        if (in) cout << infile << " is open\n";
        else cout << infile << " is not open, exiting program\n", exit(2);
        while (in.peek() == '#') in.getline(line, 1000000);
        cout << "\n\n";

        while (!in.eof()) {
            in >> temp >> locHUM >> temp >> ref >> alt >> temp >> temp >> temp >> temp;
            in.getline(line, 1000000);

            lcn2 = locHUM / 1000000;
            if (prev2 < lcn2) cout << "chomosome " << chrom << ", megabase " << prev2 << "\n";

            int lcn = locHUM / 50000;
            if (lcn > prev) { 
            }
            prev = lcn, prev2 = lcn2;

            if (temp[0] == 'G' && temp[1] == 'T') { 
                int na = 0, at = 0; // number of 
                int geno[1000][2]{}; // genotypes stored in array
                for (int j = 0; j < 934; j++) {
                    if (line[j * 4 + 1] != '.' && line[j * 4 + 3] != '.') { 
                        geno[j][0] = line[j * 4 + 1] - 48; // convert to '0' and '1'
                        geno[j][1] = line[j * 4 + 3] - 48;
                        if (j > 4) at += (geno[j][0] + geno[j][1]); 
                        na += 2; // count scored alleles
                    }
                    else geno[j][0] = -1, geno[j][1] = -1; 
                }

                // check first 5 genotypes are homozygous and scored
                bool test = true;
                for (int i = 0; i < 5; i++) if ((geno[i][0] != geno[i][1]) || (geno[i][0] == -1)) test = false;

                if (na > 1800 && test && ref[1] == '\0' && alt[1] == '\0') { 
                  int TP = 32, major = 0;
                   if (at < float(na) / 2) TP = triplet(locHUM, ref[0]); 
                   else TP = triplet(locHUM, alt[0]), major = 1;
                   int ts = TV(ref[0], alt[0]); // determine TS/TV
                   int typ = Get_type(geno, ref[0], alt[0], at); 
                   if (TP < 32) types[TP + 32 * ts][typ]++;
                  D_calculations(geno, ref[0], alt[0], TP, ts, major);                  
                }
            }
        }
    }
    Generate_output(out);
    return 0;
}

void D_calculations(int geno[1000][2], char rf, char at, int trip, int tstv, int major)
{
    float freq[54][2]{};
    float FQ[7][2]{};

    for (int i = 5; i < 934; i++) { // determine population (freq[pop][]) and regional (FQ[region][]) alternate allele counts [0] = aa, [1] = all alleles
        if (geno[i][0] == 1) freq[pops[i][0]][0]++, FQ[pops[i][1]][0]++;
        if (geno[i][1] == 1) freq[pops[i][0]][0]++, FQ[pops[i][1]][0]++;
        if (geno[i][0] != -1) freq[pops[i][0]][1] += 2, FQ[pops[i][1]][1] += 2; // only count scored alleles
    }
    trip += tstv * 32;

    float fP1, fP2;
    if (t1 == 0) fP1 = freq[P1][0] / freq[P1][1];
    else fP1 = FQ[P1][0] / FQ[P1][1];
    if (t2 == 0) fP2 = freq[P2][0] / freq[P2][1];
    else fP2 = FQ[P2][0] / FQ[P2][1];

    int typ;
    if (geno[0][0] == 0) { // chimpanzee = reference
        typ = -1;
        if (fP1 == 0 && fP2 > 0) typ = 0; // Africans fixed for A
        if (fP1 > 0 && fP2 == 0) typ = 1; // Europe fixed for A
        if (fP1 == 1 && fP2 < 1) typ = 2; // Africans fixed for B
        if (fP1 < 1 && fP2 == 1) typ = 3; // Europe fixed for B
        if (fP1 > 0 && fP1 < 1 && fP2 > 0 && fP2 < 1) typ = 4; // both poly

        for (int arch = 2; arch < 5; arch++) {
            if (geno[arch][0] == 1 && typ > -1) {
                ABBA[trip][typ + ((arch - 2) + 3 * major) * 6][0] += fP1 * (1 - fP2);
                ABBA[trip][typ + ((arch - 2) + 3 * major) * 6][1] += fP2 * (1 - fP1);
                ABBA[trip][5 + ((arch - 2) + 3 * major) * 6][0] += fP1 * (1 - fP2);
                ABBA[trip][5 + ((arch - 2) + 3 * major) * 6][1] += fP2 * (1 - fP1);
            }
        }
    }

    if (geno[0][0] == 1) { // chimpanzee = alternate
        major = 1 - major;
        typ = -1;
        if (fP1 == 0 && fP2 > 0) typ = 2; // Africans fixed for B
        if (fP1 > 0 && fP2 == 0) typ = 3; // Europe fixed for B
        if (fP1 == 1 && fP2 < 1) typ = 0; // Africans fixed for A
        if (fP1 < 1 && fP2 == 1) typ = 1; // Europe fixed for A
        if (fP1 > 0 && fP1 < 1 && fP2 > 0 && fP2 < 1) typ = 4; // both poly

        for (int arch = 2; arch < 5; arch++) {
            if (geno[arch][0] == 0 && typ > -1) {
                ABBA[trip][typ + ((arch - 2) + 3 * major) * 6][1] += fP1 * (1 - fP2);
                ABBA[trip][typ + ((arch - 2) + 3 * major) * 6][0] += fP2 * (1 - fP1);
                ABBA[trip][5 + ((arch - 2) + 3 * major) * 6][1] += fP1 * (1 - fP2);
                ABBA[trip][5 + ((arch - 2) + 3 * major) * 6][0] += fP2 * (1 - fP1);
            }
        }
    }
}

int TV(char A, char B)
{
    if ((A == 'A' && B == 'G') || (A == 'G' && B == 'A')) return 1;
    if ((A == 'C' && B == 'T') || (A == 'T' && B == 'C')) return 1;
    return 0;
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

int Get_type(int geno[1000][2], char rf, char at, int alternate) // calculates various path lengths and ABBA and BABA counts
{
    if (alternate == 0) {
        if (geno[0][0] == 0 && geno[2][0] == 1 && geno[3][0] == 0 && geno[4][0] == 0) return 0; // Vindija different
        if (geno[0][0] == 0 && geno[2][0] == 0 && geno[3][0] == 1 && geno[4][0] == 0) return 1; // Altai different
        if (geno[0][0] == 0 && geno[2][0] == 0 && geno[3][0] == 0 && geno[4][0] == 1) return 2; // Denisovan different
        if (geno[0][0] == 1 && geno[2][0] == 1 && geno[3][0] == 1 && geno[4][0] == 1) return 3; // humans different
        if (geno[0][0] == 0 && geno[2][0] == 1 && geno[3][0] == 1 && geno[4][0] == 0) return 4; // Neanderthals different
    }
    else {
        if (geno[0][0] != geno[2][0] && geno[0][0] == geno[3][0] && geno[0][0] == geno[4][0]) return 5; // Vindija different, humans poly
        if (geno[0][0] != geno[3][0] && geno[0][0] == geno[2][0] && geno[0][0] == geno[4][0]) return 6; // Altai different, humans poly
        if (geno[0][0] != geno[4][0] && geno[0][0] == geno[2][0] && geno[0][0] == geno[3][0]) return 7; // Denisovan different
        if (geno[0][0] == geno[4][0] && geno[2][0] == geno[3][0] && geno[0][0] != geno[2][0]) return 8; // Neanderthal different, humans poly
        if (geno[0][0] == geno[2][0] && geno[2][0] == geno[3][0] && geno[3][0] == geno[4][0]) return 9; // just humans poly
        if (geno[0][0] != geno[2][0] && geno[2][0] == geno[3][0] && geno[3][0] == geno[4][0]) return 10; // chimp different, humans poly
    }
    return 11;
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

void LoadPopsHGDP(std::wstring SF) // read in population and region codes for each of 2504 individuals
{
    int N_ind[60]{};
    int N_reg[60]{};
    for (int i = 0; i < 100; i++) popnams[i][0] = 'w';
    std::string infile(SF.begin(), SF.end());
    infile += "inpopsHGDP.txt";
    std::ifstream in1(infile);

    char temp[100], nam[100];
    for (int i = 0; i < 929; i++) {
        in1 >> temp >> nam >> temp >> temp >> pops[i + 5][0] >> pops[i + 5][1];
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
        cout << "end chromosome\t\t\t" << nd << "\t\t" << "\t<f>\n\n";
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
        if (k != 'a' && k != 'b' && k != 'c' && k != 'd' && k != 'e' && k != 'f' && k != 'Y') {
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
    out << "Overall counts of ABBAs and BABAs according to site class\n\n";
    out << "P1 = ";
    if (t1 == 0) out << popnams[P1] << "\n";
    else out << popnams[P2] << "\n";
    out << "P2 = ";
    if (t2 == 0) out << popnams[P2] << "\n";
    else out << popnams[P2] << "\n";
    out << "chromosomes " << st << " to " << nd << "\n\n";
    
    // add headers output columns
    out << "Count of sites by class\n";
    out << "triplet" << "\t" << "Vin dif" << "\t" << "Altai dif" << "\t" << "Den dif" << "\t" << "humans dif" << "\t" << "Neands dif";
    out << "\t" << "Vin dif HP" << "\t" << "Altai dif HP" << "\t" << "Den dif HP" << "\t" << "Neas dif HP" << "\t" << "HP" << "\t" << "Cmp dif HP" << "\n";

    for (int i = 0; i < 64; i++) {
        if (i < 32) out << triptype[i] << "_TV";
        else out << triptype[i - 32] << "_TS";
        for (int j = 0; j < 11; j++) out << "\t" << types[i][j];
        out << "\n";
    }

    out << "\n\n";
    out << popnams[P1] << "\t" << popnams[P2] << "\n";
    out << "\t" << "major = chimpanzee A" << "\t\t\t\t\t\t\t\t\t\t\t\t" << "major = chimpanzee A" << "\t\t\t\t\t\t\t\t\t\t\t\t" << "major = chimpanzee A" << "\t\t\t\t\t\t\t\t\t\t\t\t";
    out << "major = archaic B" << "\t\t\t\t\t\t\t\t\t\t\t\t" << "major = archaic B" << "\t\t\t\t\t\t\t\t\t\t\t\t" << "major = archaic B" << "\n";
    out << "\t" << popnams[70] << "\t\t\t\t\t\t\t\t\t\t\t\t" << popnams[71] << "\t\t\t\t\t\t\t\t\t\t\t\t" << popnams[72] << "\t\t\t\t\t\t\t\t\t\t\t\t";
    out << popnams[70] << "\t\t\t\t\t\t\t\t\t\t\t\t" << popnams[71] << "\t\t\t\t\t\t\t\t\t\t\t\t" << popnams[72] << "\n";
    out << "\t";
    for (int i = 0; i < 6; i++) out << "A/P" << "\t\t" << "P/A" << "\t\t" << "B/P" << "\t\t" << "P/B" << "\t\t" << "P/P" << "\t\t" << "ALL" << "\t\t";
    out << "\n\t";
    for (int i = 0; i < 36; i++) out << "BABA" << "\t" << "ABBA" << "\t";
    out << "\n";

    for (int i = 0; i < 64; i++) {
        if (i < 32) out << triptype[i] << "_TV";
        else out << triptype[i - 32] << "_TS";
        for (int j = 0; j < 36; j++) out << "\t" << ABBA[i][j][0] << "\t" << ABBA[i][j][1];
        out << "\n";
    }
    out << "\n\n";

    out << "Table2\n";
    long counts[12][2]{};
    for (int tstv = 0; tstv < 2; tstv++) {
        for (int j = 0; j < 11; j++) {
            for (int k = 0; k < 32; k++) counts[j][tstv] += types[k + tstv * 32][j];
        }
    }
    out << "AAABA" << "\t" << counts[2][0] << "\t" << counts[2][1] << "\t" << counts[2][1] / counts[2][0] << "\n";
    out << "AABAA" << "\t" << counts[1][0] << "\t" << counts[1][1] << "\t" << counts[1][1] / counts[1][0] << "\n";
    out << "ABAAA" << "\t" << counts[0][0] << "\t" << counts[0][1] << "\t" << counts[0][1] / counts[0][0] << "\n";
    out << "BAAAA" << "\t" << counts[3][0] << "\t" << counts[3][1] << "\t" << counts[3][1] / counts[3][0] << "\n";
    out << "ABBAA" << "\t" << counts[4][0] << "\t" << counts[4][1] << "\t" << counts[4][1] / counts[4][0] << "\n";

    out << "PAAAA" << "\t" << counts[9][0] << "\t" << counts[9][1] << "\t" << counts[9][1] / counts[9][0] << "\n";
    out << "PAABA" << "\t" << counts[7][0] << "\t" << counts[7][1] << "\t" << counts[7][1] / counts[7][0] << "\n";
    out << "PABAA" << "\t" << counts[6][0] << "\t" << counts[6][1] << "\t" << counts[6][1] / counts[6][0] << "\n";
    out << "PBAAA" << "\t" << counts[5][0] << "\t" << counts[5][1] << "\t" << counts[5][1] / counts[5][0] << "\n";
    out << "PBBBA" << "\t" << counts[10][0] << "\t" << counts[10][1] << "\t" << counts[10][1] / counts[10][0] << "\n";
    out << "PBBAA" << "\t" << counts[8][0] << "\t" << counts[8][1] << "\t" << counts[8][1] / counts[8][0] << "\n";
    out << "\n\n";

    out << "Figure 7\n";
    out << "\t" << "AABAA" << "\t\t" << "PAAAA" << "\t\t" << "PBBBA" << "\n";
    out << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\n";
    for (int k = 0; k < 32; k++) {
        out << triptype[k];
        out << "\t" << float(types[k][1]) / counts[1][0] << "\t" << float(types[k + 32][1]) / counts[1][1];
        out << "\t" << float(types[k][9]) / counts[9][0] << "\t" << float(types[k + 32][9]) / counts[9][1];
        out << "\t" << float(types[k][10]) / counts[10][0] << "\t" << float(types[k + 32][10]) / counts[10][1] << "\n";
    }
    out << "\n\n";

    out << "Figure 8" << "\n";
    out << "Triplet" << "\t" << "propn ABBA" << "\t" << "propn BABA" << "\t" "difference" << "\t" << "propn PBBBA" << "\t" << "propn PAAAA" << "\t" << "difference" << "\n";
    float AB = 0, BA = 0;
    for (int k = 0; k < 32; k++) {
        AB += (ABBA[k][11][1] + ABBA[k][29][1]);
        BA += (ABBA[k][11][0] + ABBA[k][29][0]);
    }
    for (int k = 0; k < 32; k++) {
        out << triptype[k];
        out << "\t" << float(ABBA[k][11][1] + ABBA[k][29][1]) / AB << "\t" << float(ABBA[k][11][0] + ABBA[k][29][0]) / BA;
        out << "\t" << float(ABBA[k][11][1] + ABBA[k][29][1]) / AB - float(ABBA[k][11][0] + ABBA[k][29][0]) / BA;
        out << "\t" << float(types[k][10]) / counts[10][0] << "\t" << float(types[k][9]) / counts[9][0] << "\t" << float(types[k][10]) / counts[10][0] - float(types[k][9]) / counts[9][0] << "\n";
    }
    out << "\n\n";

    out << "Figure 9" << "\n";
    out << "\t" << "Tranversions" << "\t\t" << "Transitions\n";
    out << "Triplet" << "\t" << "D major A" << "\t" << "D major B" << "\t" << "D major A" << "\t" << "D major B" << "\n";
    for (int k = 0; k < 32; k++) {
        out << triptype[k];
        out << "\t" << float(ABBA[k][11][1] - ABBA[k][11][0]) / (ABBA[k][11][1] + ABBA[k][11][0]) << "\t" << float(ABBA[k][29][1] - ABBA[k][29][0]) / (ABBA[k][29][1] + ABBA[k][29][0]);
        out << "\t" << float(ABBA[k + 32][11][1] - ABBA[k + 32][11][0]) / (ABBA[k + 32][11][1] + ABBA[k + 32][11][0]);
        out << "\t" << float(ABBA[k + 32][29][1] - ABBA[k + 32][29][0]) / (ABBA[k + 32][29][1] + ABBA[k + 32][29][0]) << "\n";
    }
    out << "\n\n";

    out << "Table 3" << "\n";
    out << "\t" << "Altai" << "\t\t" << "Vindija" << "\t\t" << "Denisovan\n";
    out << "P1/P2" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\n";
    float Dstar[6][3][2]{};
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 6; k++) {
            for (int arch = 0; arch < 3; arch++) {
                Dstar[k][arch][0] += (ABBA[j][arch * 6 + k][1] - ABBA[j][arch * 6 + k][0] + ABBA[j][arch * 6 + k + 18][1] - ABBA[j][arch * 6 + k + 18][0]);
                Dstar[k][arch][1] += (ABBA[j + 32][arch * 6 + k][1] - ABBA[j + 32][arch * 6 + k][0] + ABBA[j + 32][arch * 6 + k + 18][1] - ABBA[j + 32][arch * 6 + k + 18][0]);
            }
        }
    }
    float T0 = Dstar[5][0][0] + Dstar[5][0][1], T1 = Dstar[5][1][0] + Dstar[5][1][1], T2 = Dstar[5][2][0] + Dstar[5][2][1];
    out << "A/P" << "\t" << 100 * Dstar[0][0][0] / T0 << "\t" << 100 * Dstar[0][0][1] / T0 << "\t" << 100 * Dstar[0][1][0] / T1 << "\t" << 100 * Dstar[0][1][1] / T1 << "\t" << 100 * Dstar[0][2][0] / T2 << "\t" << 100 * Dstar[0][2][1] / T2 << "\n";
    out << "P/A" << "\t" << 100 * Dstar[1][0][0] / T0 << "\t" << 100 * Dstar[1][0][1] / T0 << "\t" << 100 * Dstar[1][1][0] / T1 << "\t" << 100 * Dstar[1][1][1] / T1 << "\t" << 100 * Dstar[1][2][0] / T2 << "\t" << 100 * Dstar[1][2][1] / T2 << "\n";
    out << "B/P" << "\t" << 100 * Dstar[2][0][0] / T0 << "\t" << 100 * Dstar[2][0][1] / T0 << "\t" << 100 * Dstar[2][1][0] / T1 << "\t" << 100 * Dstar[2][1][1] / T1 << "\t" << 100 * Dstar[2][2][0] / T2 << "\t" << 100 * Dstar[2][2][1] / T2 << "\n";
    out << "P/B" << "\t" << 100 * Dstar[3][0][0] / T0 << "\t" << 100 * Dstar[3][0][1] / T0 << "\t" << 100 * Dstar[3][1][0] / T1 << "\t" << 100 * Dstar[3][1][1] / T1 << "\t" << 100 * Dstar[3][2][0] / T2 << "\t" << 100 * Dstar[3][2][1] / T2 << "\n";
    out << "P/P" << "\t" << 100 * Dstar[4][0][0] / T0 << "\t" << 100 * Dstar[4][0][1] / T0 << "\t" << 100 * Dstar[4][1][0] / T1 << "\t" << 100 * Dstar[4][1][1] / T1 << "\t" << 100 * Dstar[4][2][0] / T2 << "\t" << 100 * Dstar[4][2][1] / T2 << "\n";
    out << "ALL" << "\t" << 100 * Dstar[5][0][0] / T0 << "\t" << 100 * Dstar[5][0][1] / T0 << "\t" << 100 * Dstar[5][1][0] / T1 << "\t" << 100 * Dstar[5][1][1] / T1 << "\t" << 100 * Dstar[5][2][0] / T2 << "\t" << 100 * Dstar[5][2][1] / T2 << "\n";
    out << "\n\n";

    out << "Figure 9" << "\n";
    out << "\t" << "A/P" << "\t\t" << "P/A" << "\t\t" << "B/P" << "\t\t" << "P/B" << "\t\t" << "P/P" << "\t\t" << "ALL" << "\t\t" << "ALL\n";
    out << "Triplet" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TV" << "\t" << "TS" << "\t" << "TVTS\n";
    for (int j = 0; j < 32; j++) {
        out << triptype[j];
        for (int k = 0; k < 6; k++) {
            out << "\t" << ABBA[j][6 + k][1] - ABBA[j][6 + k][0] + ABBA[j][6 + k + 18][1] - ABBA[j][6 + k + 18][0] << "\t" << ABBA[j + 32][6 + k][1] - ABBA[j + 32][6 + k][0] + ABBA[j + 32][6 + k + 18][1] - ABBA[j + 32][6 + k + 18][0];
        }
        out << "\t" << ABBA[j][11][1] - ABBA[j][11][0] + ABBA[j][11 + 18][1] - ABBA[j][11 + 18][0] + ABBA[j + 32][11][1] - ABBA[j + 32][11][0] + ABBA[j + 32][11 + 18][1] - ABBA[j + 32][11 + 18][0];
        out << "\n";
    }
    out.close();
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
