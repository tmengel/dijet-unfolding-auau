// Rebuild jer_smear_functions.root, replacing the hard-coded p+p base smearing
// constants with the updated PPG-08 p+p IAN values (Table 5, R = 0.3) and
// optionally relaxing the UE quadrature clamp.
//
// The stored functions are
//     sigma_extra(pT) = sqrt( BASE^2 + max(0, JER(b_var)^2 - JER(b_MC)^2) )
//     JER(b)          = sqrt( c^2 + a^2/pT + b^2/pT^2 )
// with parameters  p0,p3 = c   p1,p4 = a   p2 = b_var   p5 = b_MC.
// BASE is the only quantity hard-coded in the formula string, so the
// centrality-dependent UE terms carry over untouched from the input file.
//
// The max(0, ...) clamp fires whenever b_var <= b_MC. For 0-10% that is true for
// nominal (6.1597) and negative (5.9223) but not positive (6.397), which leaves
// the UE part of the JER systematic one-sided in the most central interval.
// unclamp = true instead floors the total at zero, letting a UE deficit reduce
// the extra smearing. Flooring the total (not the difference) is required: for
// negJER at 7 GeV the difference alone would drive the argument of the sqrt
// negative.
//
// Usage:
//   root -l -b -q 'jer/remake_smear_functions.C()'
//   root -l -b -q 'jer/remake_smear_functions.C("in.root","out.root",true)'

const double kBaseNominal  = 0.107;  // p+p IAN Table 5, R = 0.3 nominal smearing
const double kBaseNegative = 0.083;  // nominal - 2.4%
const double kBasePositive = 0.131;  // nominal + 2.4%

void remake_smear_functions(
  const char *infile  = "jer/jer_smear_functions_prelim_8pct_base.root",
  const char *outfile = "jer/jer_smear_functions.root",
  bool unclamp = false)
{
  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie())
    {
      std::cerr << "Cannot open " << infile << std::endl;
      return;
    }

  const char *var = "sqrt([0]*[0]+([1]*[1]/(x))+([2]*[2]/((x)*(x))))";
  const char *ref = "sqrt([3]*[3]+([4]*[4]/(x))+([5]*[5]/((x)*(x))))";

  TFile *fout = TFile::Open(outfile, "RECREATE");
  int nwritten = 0;

  TIter next(fin->GetListOfKeys());
  TKey *key = nullptr;
  while ((key = (TKey*) next()))
    {
      TF1 *fold = (TF1*) key->ReadObj();
      if (!fold->InheritsFrom("TF1")) continue;

      TString name = fold->GetName();
      double base = name.BeginsWith("f_nominal")  ? kBaseNominal
                  : name.BeginsWith("f_negative") ? kBaseNegative
                  : name.BeginsWith("f_positive") ? kBasePositive : -1;
      if (base < 0)
        {
          std::cerr << "Unrecognized function " << name << ", skipping" << std::endl;
          continue;
        }

      // difference of squares, either clamped on its own or folded into the total
      TString diff = Form("TMath::Power((%s),2)-TMath::Power((%s),2)", var, ref);
      TString arg  = unclamp
        ? Form("%.6f*%.6f+%s>0?%.6f*%.6f+%s:0.0", base, base, diff.Data(), base, base, diff.Data())
        : Form("%.6f*%.6f+(%s>0?%s:0.0)", base, base, diff.Data(), diff.Data());

      double xlo = 0, xhi = 0;
      fold->GetRange(xlo, xhi);
      TF1 *fnew = new TF1(name, Form("sqrt(%s)", arg.Data()), xlo, xhi);
      for (int p = 0; p < fold->GetNpar(); p++) fnew->SetParameter(p, fold->GetParameter(p));

      fout->cd();
      fnew->Write(name);
      nwritten++;

      printf("%-16s base %.3f  sigma(10,30,60) = %.4f %.4f %.4f   (was %.4f %.4f %.4f)\n",
             name.Data(), base,
             fnew->Eval(10), fnew->Eval(30), fnew->Eval(60),
             fold->Eval(10), fold->Eval(30), fold->Eval(60));
    }

  fout->Close();
  fin->Close();
  printf("\nWrote %d functions to %s%s\n", nwritten, outfile,
         unclamp ? "  (UE clamp relaxed: total floored at zero)" : "");
}
