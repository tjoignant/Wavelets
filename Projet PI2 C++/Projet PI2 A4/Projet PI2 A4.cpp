// Projet PI2 A4.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "math.h"
using namespace std;

double scalingFunction_Haar(double j, double k, double t)
{
    double coef1 = pow(2, j/2);
    double coef2 = pow(2, j);
    double test = coef2 * t - k;
    if (test <= 1 && test > 0)
    {
        return coef1;
    }
    else
    {
        return 0;
    }
}

double scalingCoef_Haar(vector<double> z, double j, double k)
{
    double size = z.size();
    double coef = 0;
    for (int t = 0; t <= min((1+k)/pow(2,j),size - 1); t++)
    {
        coef += z[t] * scalingFunction_Haar(j, k, t);
    }
    return coef;
}

double waveletsHaar(double j, double k, double t)
{
    double coef1 = std::pow(2, j / 2);
    double coef2 = pow(2, j);
    double test;
    test = coef2 * t - k;
    if (test >= 0 && test < 1)
    {
        if (test < 1.0 / 2)
        {
            return coef1;
        }
        else
        {
            return -coef1;
        }
    }
    else
    {
        return 0;
    }
}

double min(double m, double n)
{
    if (m <= n) { return m; }
    else { return n; }

}

double max(double m, double n)
{
    if (m >= n) { return m; }
    else { return n; }

}

double waveletsTransfo_Haar(vector<double> z, double j, double k, int tmax)
{
    double result = 0;
    for (int t = max((k)/pow(2,j),0); t <= min((1+k)/pow(2,j),tmax); t++)
    {
        result += z[t] * waveletsHaar(j, k, t);
    }
    return result;
}

//int main()
//{
//    vector<double> monVecteur;
//    ifstream monFlux("C:/Users/gaelm/source/repos/Projet PI2 A4/EURNOK fermeture.csv");
//    string ligne;
//    double nb;
//    char* p;
//    int i = 0;
//    if (monFlux)
//    {
//        while (getline(monFlux, ligne))
//        {
//            monVecteur.push_back(strtod(ligne.c_str(), &p));
//            i = i + 1;
//        }
//    }
//    else
//    {
//        cout << "ERREUR.";
//    }
//
//    //cout << monVecteur[0] << endl << monVecteur[5] << endl;
//    //cout << "Cas j=-6 et k=1:" << std::endl << waveletsTransfo_Haar(monVecteur, -6, 1, monVecteur.size());
//
//    monFlux.close();
//
//    double sum;
//    int size = monVecteur.size();
//    string const nomFichier("C:/Users/gaelm/source/repos/Projet PI2 A4/Dj EURNOK.csv");
//    ofstream Dj(nomFichier.c_str());
//    if (Dj)
//    {
//        for (int t = 1; t <= size; t++)
//        {
//            for (int j = -5; j <= 5; j++)
//            {
//                sum = 0;
//                for (int k = -size; k <= size; k++)
//                {
//                    sum += waveletsTransfo_Haar(monVecteur, j, k, size - 1) * waveletsHaar(j, k, t);
//                }
//                if (j == 5)
//                {
//                    Dj << sum << endl;
//                }
//                else
//                {
//                    Dj << sum << ";";
//                }
//            }
//        }
//    }
//    else
//    {
//        cout << "ERREUR";
//    }
//    Dj.close();
//
//}

//Test de l'approximation à l'aide de la fonction de Haar

//int main()    
//{
//    vector<double> monVecteur;
//    ifstream monFlux("C:/Users/gaelm/source/repos/Projet PI2 A4/EURNOK fermeture.csv");
//    string ligne;
//    double nb;
//    char* p;
//    int i = 0;
//    if (monFlux)
//    {
//        while (getline(monFlux, ligne))
//        {
//            monVecteur.push_back(strtod(ligne.c_str(), &p));
//            i = i + 1;
//        }
//    }
//    else
//    {
//        cout << "ERREUR.";
//    }
//
//
//    monFlux.close();
//
//    double sum;
//    int size = monVecteur.size();
//    string const nomFichier("C:/Users/gaelm/source/repos/Projet PI2 A4/Approx EURNOK.csv");
//    ofstream Approx(nomFichier.c_str());
//    int j = -2;
//    if (Approx)
//    {
//        for (int t = 1; t < size; t++)
//        {
//            sum = 0;
//            for (int k = -size; k <= size; k++)
//            {
//                sum += scalingCoef_Haar(monVecteur, j, k) * scalingFunction_Haar(j, k, t);
//            }
//            Approx << monVecteur[t] << ";" << sum << endl;
//        }
//    }
//    else
//    {
//        cout << "ERREUR";
//    }
//    Approx.close();
//
//}

 //Test de l'approximation à l'aide de la fonction de Haar avec une méthode plus rapide

//int main()   
//{
//    vector<double> monVecteur;
//    ifstream monFlux("C:/Users/gaelm/source/repos/Projet PI2 A4/EURNOK fermeture.csv");
//    string ligne;
//    double nb;
//    char* p;
//    int i = 0;
//    if (monFlux)
//    {
//        while (getline(monFlux, ligne))
//        {
//            monVecteur.push_back(strtod(ligne.c_str(), &p));
//            i = i + 1;
//        }
//    }
//    else
//    {
//        cout << "ERREUR.";
//    }
//
//
//    monFlux.close();
//
//    int size = monVecteur.size();
//    vector<double> sum = vector<double>(size,0);
//    double function;
//    string const nomFichier("C:/Users/gaelm/source/repos/Projet PI2 A4/Approx EURNOK.csv");
//    ofstream Approx(nomFichier.c_str());
//    int j = -2;
//    if (Approx)
//    {
//        for (int k = -pow(2, j) * size; k <= pow(2, j) * size - 1; k++)
//        {
//            for (int t = 1; t < size; t++)
//            {
//                function = scalingFunction_Haar(j, k, t);
//                sum[t] += scalingCoef_Haar(monVecteur, j, k) * function;
//            }
//            cout << "Cas k = " << k << endl;
//        }
//        for (int t = 1; t < size; t++)
//        {
//            Approx << monVecteur[t] << ";" << sum[t] << endl;
//        }
//    }
//    else
//    {
//        cout << "ERREUR";
//    }
//    Approx.close();
//
//}


vector<double> coefficients_Daubechies(int moment)
{
    vector<double> coef;
    switch (moment)
    {
    case 1:
        coef.push_back(1);
        coef.push_back(1);
        return coef;
    case 2:
        coef.push_back(0.6830127);
        coef.push_back(1.1830127);
        coef.push_back(0.3169873);
        coef.push_back(-0.1830127);
        return coef;
    case 3:
        coef.push_back(0.47046721);
        coef.push_back(1.14111692);
        coef.push_back(0.650365);
        coef.push_back(-0.19093442);
        coef.push_back(-0.12083221);
        coef.push_back(0.0498175);
        return coef;
    case 4:
        coef.push_back(0.32580343);
        coef.push_back(1.01094572);
        coef.push_back(0.89220014);
        coef.push_back(-0.03957503);
        coef.push_back(-0.26450717);
        coef.push_back(0.0436163);
        coef.push_back(0.0465036);
        coef.push_back(-0.01498699);
        return coef;
    }
    coef.push_back(0);
    return coef;
}

double getPhi(double t, vector<double> t_tampon, vector<double> phi_t_tampon)
{
    double phi = 0;
    for (int i = 0; i < t_tampon.size(); i++)
    {
        if (t == t_tampon[i])
        {
            phi = phi_t_tampon[i];
            break;
        }
    }
    return phi;
}

vector<vector<double>> buildFatherWavelet(int moment, int depth)
{
    vector<vector<double>> result;
    vector<double> support;
    if (moment == 1)
    {
        support.push_back(0);
        support.push_back(1);
    }
    else
    {
        support.push_back(0);
        support.push_back(10);
    }
    int step = 1;
    vector<double> t(support[1] - support[0] + 1, 0);
    vector<double> phi_t(support[1] - support[0] + 1, 0);
    for (int i = 0; i < t.size(); i++)
    {
        t[i] = support[0] + i;
        if (t[i] == 1)
        {
            phi_t[i] = 1;
        }
    }
    vector<double> t_tampon;
    vector<double> phi_t_tampon;
    double sum;
    vector<double> coef = coefficients_Daubechies(moment);
    // Cascade algorithm
    for (int j = 2; j < depth + 1; j++)
    {
        step = step / 2;
        t_tampon = t;
        phi_t_tampon = phi_t;
        t.clear();
        t = vector<double>(2*t_tampon.size()-1, 0);
        phi_t.clear();
        phi_t = vector<double>(2 * phi_t_tampon.size() - 1, 0);
        for (int i = 0; i < t.size(); i = i + 2)
        {
            t[i] = t_tampon[i / 2];
            phi_t[i] = phi_t_tampon[i / 2];
        }
        for (int i = 1; i < t.size(); i = i + 2)
        {
            t[i] = (t_tampon[(i / 2)] + t_tampon[(i / 2) + 1]) / 2.0;
            sum = 0;
            for (int k = 0; k < 2 * moment; k++)
            {
                sum += coef[k] * getPhi(2 * t[i] - k, t_tampon, phi_t_tampon);
            }
            phi_t[i] = sum;
        }
    }
    // On retire tous les 0 lorsqu'il n'y a plus de valeur après.
    int index = 0;
    bool check = false;
    for (int i = 0; i < phi_t.size(); i++)
    {
        if (phi_t[i] == 0 && !check)
        {
            index = i;
            check = true;
        }
        if (phi_t[i] != 0 && check)
        {
            check = false;
        }
    }
    for (int i = phi_t.size() - 1; i >= index; i--)
    {
        t.pop_back();
        phi_t.pop_back();
    }
    result.push_back(t);
    result.push_back(phi_t);
    return result;
}

// Test de la scaling function de daubechies.

//int main()  
//{
//
//    int moment = 2;
//    int depth = 6;
//
//    vector<vector<double>> building = buildFatherWavelet(moment, depth);
//    vector<double> time = building[0];
//    vector<double> fatherWavelet = building[1];
//
//    double sum;
//    int size = time.size();
//    string const nomFichier("C:/Users/gaelm/source/repos/Projet PI2 A4/Approx EURNOK.csv");
//    ofstream Approx(nomFichier.c_str());
//    int j = -2;
//    if (Approx)
//    {
//        for (int t = 1; t < size; t++)
//        {
//            Approx << time[t] << ";" << fatherWavelet[t] << endl;
//        }
//    }
//    else
//    {
//        cout << "ERREUR";
//    }
//    Approx.close();
//
//}



double scalingFunction_Daubechies(double j, double k, double t, vector<vector<double>> building)
{
    double phi = 0;
    vector<double> time = building[0];
    vector<double> fatherWavelet = building[1];
    double coef1 = pow(2, j / 2);
    double coef2 = pow(2, j);
    double test = coef2 * t - k;
    if (test >= time[0] && test <= time[time.size() - 1])
    {
        for (int i = 0; i < time.size(); i++)
        {
            if (test == time[i])
            {
                phi = coef1 * fatherWavelet[i];
                break;
            }
            else if (test < time[i])
            {
                phi = coef1 * ((fatherWavelet[i - 1] - fatherWavelet[i]) * t + time[i - 1] * fatherWavelet[i] - time[i] * fatherWavelet[i - 1]) / (time[i - 1] - time[i]);  // interpolation
                break;
            }
        }
    }
    return phi;
}

double scalingCoef_Daubechies(vector<double> z, double j, double k, vector<vector<double>> building)
{
    int size = z.size();
    double coef = 0;
    for (int t = 0; t < size; t++)
    {
        coef += z[t] * scalingFunction_Daubechies(j, k, t, building);
    }
    return coef;
}

// Test de l'approximation avec Daubechies

int main() 
{
    vector<double> monVecteur;
    ifstream monFlux("C:/Users/gaelm/source/repos/Projet PI2 A4/EURNOK fermeture.csv");
    string ligne;
    double nb;
    char* p;
    int i = 0;
    if (monFlux)
    {
        while (getline(monFlux, ligne))
        {
            monVecteur.push_back(strtod(ligne.c_str(), &p));
            i = i + 1;
        }
    }
    else
    {
        cout << "ERREUR.";
    }


    monFlux.close();

    int moment = 4;
    int depth = 10;
    int j = -1;

    vector<vector<double>> building = buildFatherWavelet(moment, depth);

    int size = monVecteur.size();
    vector<double> sum = vector<double>(size, 0);

    vector<double> coef_Daubechies;
    for (int k = -size; k <= size; k++)
    {
        coef_Daubechies.push_back(scalingCoef_Daubechies(monVecteur, j, k, building));
    }

    string const nomFichier("C:/Users/gaelm/source/repos/Projet PI2 A4/Approx1111 EURNOK.csv");
    ofstream Approx(nomFichier.c_str());
    if (Approx)
    {
        for (int k = -pow(2, j) * size; k <= pow(2, j) * size - 1; k++)
        {
            for (int t = 1; t < size; t++)
            {
                sum[t] += coef_Daubechies[k + size] * scalingFunction_Daubechies(j, k, t, building);
            }
            cout << "Cas k = " << k << endl;
        }
        for (int t = 1; t < size; t++)
        {
            Approx << monVecteur[t] << ";" << sum[t] << endl;
        }
    }
    else
    {
        cout << "ERREUR";
    }
    Approx.close();

}

// Estimation fonction de densité

