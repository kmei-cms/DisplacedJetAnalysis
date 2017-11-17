namespace LHAPDF {
    void     initPDFSet(int nset, const std::string& filename, int member=0);
	int      numberPDF(int nset);
	void     usePDFMember(int nset, int member);
	double   xfx(int nset, double x, double Q, int fl);
	double   getXmin(int nset, int member);
	double   getXmax(int nset, int member);
	double   getQ2min(int nset, int member);
	double   getQ2max(int nset, int member);
	void     extrapolate(bool extrapolate=true);
}
