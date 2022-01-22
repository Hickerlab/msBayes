%define debug_package  %{nil}

%define name msbayes
%define version 20220121
%define release 1

%define prefix   /usr/local
%define builddir $RPM_BUILD_DIR/msBayes-master

Summary: A program for testing comparative phylogeographic histories
Name: %{name}
Version: %{version}
Release: %{release}
Group: Applications/Scientific
License: GPL
Packager: Naoki Takebayashi <ntakebayashi@alaska..edu>
URL: https://github.com/Hickerlab/msBayes
Source0: msBayes-master.zip
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
Requires: R-core
BuildRequires: gsl-devel

%description

msBayes is a pipeline for testing comparative phylogeographic
histories using hierarchical ABC, developped by Michael Hickerson.  It
can be used to test for simultaneous divergence (TSD) of multiple
taxon pairs.

For the full description, see:

Overcast, I., J. C. Bagley, & M. J. Hickerson. 2017. Strategies for
  improving approximate Bayesian computation tests for synchronous
  diversification. BMC Evolutionary Biology, 17: 1-11.

Huang, W., N. Takebayashi, Y. Qi, M.J. Hickerson. 2011. MTML-msBayes:
  Approximate Bayesian compartive phylogeographic inference from
  multiple taxa and multiple loci with rate heterogeneity.  BMC
  Bioinformatics 12:1-14

Hickerson, M.J., E. Stahl, and H.A. Lessios. 2006. Test for
simultaneous divergence using approximate Bayesian computation
(ABC). Evolution 60: 2435-2453

%prep
%setup -n msBayes-master

%build
cd src; make 

%install
function CheckBuildRoot() {
    # do a few sanity checks on the BuildRoot
    # to make sure we don't damage a system
    case "${RPM_BUILD_ROOT}" in
        ''|' '|/|/bin|/boot|/dev|/etc|/home|/lib|/mnt|/root|/sbin|/tmp|/usr|/var)
            echo "Yikes!  Don't use '${RPM_BUILD_ROOT}' for a BuildRoot!"
            echo "The BuildRoot gets deleted when this package is rebuilt;"
            echo "something like '/tmp/build-blah' is a better choice."
            return 1
            ;;
        *)  return 0
            ;;
    esac
}

function CleanBuildRoot() {
    if CheckBuildRoot; then
	rm -rf "${RPM_BUILD_ROOT}"
    else
        exit 1
    fi
}

CleanBuildRoot

cd src; make PREFIX=${RPM_BUILD_ROOT}/%{prefix} install

%clean
rm -r $RPM_BUILD_ROOT
rm -r %{builddir}

%files
%defattr(-,root,root)
%doc documents/* INSTALL LICENSE README.md VERSION scripts examples
%{prefix}/bin/*
%{prefix}/lib/msbayes

%changelog
* Fri Jan 21 2022 Naoki Takebayashi <ntakebayashi@alaska.edu> [20220121-1]
- version update
- removed LDFLAGS='-static' since gsl-static is gone from Fedora

* Wed Mar  5 2014 Naoki Takebayashi <ntakebayashi@alaska.edu> [20140305-1]
- minor version update to deal with newer VGAM

* Tue Jun 11 2013 Naoki Takebayashi <ntakebayashi@alaska.edu> [20130611-1]
- version update
- minor bug fix

* Wed May 10 2012 Naoki Takebayashi <ntakebayashi@alaska.edu> [20120510-1]
- version update
- Fixed a bug that msbayes.pl get stuck without completing all runs.

* Wed Feb 22 2012 Naoki Takebayashi <ntakebayashi@alaska.edu> [20120222-1]
- version update

* Wed May 19 2010 Naoki Takebayashi <ntakebayashi@alaska.edu> [20100519-1]
- version update

* Wed May 05 2010 Naoki Takebayashi <ffnt@uaf.edu> [20100506-1]
- bug related to acceptRej.pl is fixed

* Wed May 05 2010 Naoki Takebayashi <ffnt@uaf.edu> [20100505-1]
- version update
- multi-locus version

* Thu Nov 06 2008 Naoki Takebayashi <ffnt@uaf.edu> [20081106-1]
- version update

* Fri Jun 13 2008 Naoki Takebayashi <ffnt@uaf.edu> [20080613-1]
- version update

* Thu May 15 2008 Naoki Takebayashi <ffnt@uaf.edu> [20080515-1]
- version update

* Mon Nov  6 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061106-1]
- version update

* Fri Oct 13 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061003-2]
- compiling with static library.

* Fri Oct 13 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061005-1]
- new version with correct example files

* Tue Oct 03 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061003-1]
- first release
