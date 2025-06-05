// This file contains functions related to membranes.

/*----------------------------- Calculate force ------------------------------*/
void CalcMembraneSpringForces(void) {
  int m, n, k, l, ind, side, kind, locInd, nBd, *chkBondForce;
  int nextInd, locNextInd, chkErr;
  int sideAct;
  double stf, len, eqLen, dr[NDIM], f, *membSprStf, conc;
  MembMyoCyto *pM; 
  Force *pS;

  nBd = nChMb - nMbAct;
  chkBondForce = allIntL;
  memset(chkBondForce, -1, sizeof(int) * nMbMe * nBd);
  // Spring forces between membrane points
  FOR_MBME(m) {
	kind = (ISNUC(memb.idx[m])) ? 1 : 0;
	for(n = 0; n < nBd; n++) {
		ind = P2A(memb.ch,m,n,nChMb);
		CONT(ind < 0);
		locInd = iMb[ind];
		CONT(P2A(chkBondForce,m,n,nBd) != -1);
		len = CalcVecDist(dr, &P2(memb.r,m,0), &P2(memb.r,locInd,0), 0);
		// Find the equilibrium length of a chain
		if (dimMbNuc == 2) {
			eqLen = (kind == 0) ? memb.len : memb.nucLen;
		}
		else {
			eqLen = P2A(memb.eqLen,m,n,nBd);
		}
		// Normalize the chain length using the equilibrium one
		len /= eqLen;
		pS = (kind == 0) ? &memb.spr : &memb.nucSpr;
		if (len < pS->lo || len > pS->hi) {
			f = SPRING(pS->stf, len, ((len < pS->lo) ? pS->lo : pS->hi));
			chkErr = CheckLargeForce(f, 14);
			if (chkErr != 1) {
				RecordErrorSpringForce(memb.id[m], ind, f, len, 3);
			}
			AddSpringForce(f, 1., dr, &P2(memb.f,m,0),
					&P2(memb.f,locInd,0));
    		recMb.sprF[m] += REVSIGN(f);
			if (locInd < nMbMe) {
 				recMb.sprF[locInd] += REVSIGN(f);
 			}
		}
		CONT(!(locInd < nMbMe));
		side = FindElementArray(&P2A(memb.ch,locInd,0,nChMb),
				nBd, memb.id[m], 0, 1);
		P2A(chkBondForce,locInd,side,nBd) = 1;
	}
  }

  // Spring forces between membrane points and the bound actins
  FOR_MBMECP(m) {
	pS = ISNUC(memb.idx[m]) ? &memb.nucSpr2 : &memb.spr2;
	for(n = 0; n < nMbAct; n++) {
		ind = P2A(memb.ch,m,n + nBd,nChMb);
		CONT(ind < 0);
		locInd = iAct[ind];
		// If both of the membrane point and actin belong to other subdomains,
		// the calculation doesn't need to be done here.
		CONT(m >= nMbMe && (locInd < 0 || locInd >= nActMe));
		CalcVec(dr, &P2(memb.r,m,0), &P2(act.r,locInd,0)); 
		if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }


		sideAct = FindElementArray(&P2A(act.mbReb,locInd,0,nRebMb), 
				nRebMb, memb.id[m], 0, 1);
		len = V3LEN(dr);
		f = SPRING(pS->stf, len, P2A(act.mbRebEq,locInd,sideAct,nRebMb));
		chkErr = CheckLargeForce(f, 15);
		if (chkErr != 1) {
			RecordErrorSpringForce(memb.id[m], ind, f, len, 4);
		}
		AddSpringForce(f, len, dr, &P2(memb.f,m,0), &P2(act.f,locInd,0));
	}
  }
}

void CalcMembraneBendForces(void) {
  int n, k, ind1, ind2, locInd1, locInd2, *pArr;
  double dr1[NDIM], dr2[NDIM], f1[NDIM], f2[NDIM], eqAng, stf, f;

  FOR_MBMECP(n) {
	pArr = &P2A(memb.ch,n,0,nChMb);
	for(k = 0; k < 6; k++) {
		ind1 = pArr[(k + 5) % 6];
		ind2 = pArr[(k + 1) % 6];
		CONT(ind1 < 0 || ind2 < 0);
		locInd1 = iMb[ind1];
		locInd2 = iMb[ind2];
		CONT(n >= nMbMe && !((locInd1 < nMbMe || locInd2 < nMbMe) 
				&& (locInd1 > -1 && locInd2 > -1)));
		eqAng = PI / 3.;
		stf = memb.bend.stf;
		CalcVec(dr1, &P2(memb.r,n,0), &P2(memb.r,locInd1,0));	
		CalcVec(dr2, &P2(memb.r,locInd2,0), &P2(memb.r,n,0));
		f = CalcBendForce(dr1, dr2, &P2(memb.f,locInd1,0), &P2(memb.f,n,0), 
				&P2(memb.f,locInd2,0), stf, eqAng, f1, f2);
	}
  }
}

void Calc3dMembraneBendForces(void) {
  int n, k, l, loc, CS, *pArr, ind[4], listInd[2];
  double dr[6][NDIM], crs[4][NDIM], mag[3], e[NDIM], fi[4][NDIM], tempDbl[NDIM];
  double f, c, b11, b12, b22;
  ListInt list;

  MALLOC(list.l, int, memb.unit.c * 2 * 2);
  list.c = 0;
  for(n = 0; n < memb.unit.c; n++) {
	pArr = &P2A(memb.unit.l,n,0,3);
	for(k = 0; k < 3; k++) {
		CalcVec(dr[k], &P2(memb.r,iMb[pArr[(k + 1) % NDIM]],0), 
				&P2(memb.r,iMb[pArr[k]],0));
	}
	for(k = 0; k < 3; k++) {
		ind[0] = pArr[k];
		ind[1] = pArr[(k + 1) % NDIM];
		ind[2] = pArr[(k + 2) % NDIM];
		loc = FindElementArray(&P2A(memb.ch,iMb[ind[1]],0,nChMb), 
				nChMb - nMbAct, ind[0], 0, 1);
		loc = (loc + 4) % 6;

		ind[3] = P2A(memb.ch,iMb[ind[1]],loc,nChMb);

		if (ind[0] < ind[3]) { V2SET(listInd, ind[0], ind[3]); }
		else { V2SET(listInd, ind[3], ind[0]); }
		CS = Find2ElementArray(list.l, list.c, listInd[0], listInd[1], 0, 2);
		CONT(CS != -1);
		VS3COPY(dr[3], dr[(k + 2) % NDIM], -1.);
		V3CROSS(crs[0], dr[k], dr[3]);
		mag[0] = V3LEN_SQ(crs[0]);		

		CalcVec(dr[4], &P2(memb.r,iMb[ind[2]],0), &P2(memb.r,iMb[ind[3]],0));
		CalcVec(dr[5], &P2(memb.r,iMb[ind[1]],0), &P2(memb.r,iMb[ind[3]],0));
		V3CROSS(crs[1], dr[4], dr[5]);
		mag[1] = V3LEN_SQ(crs[1]);		
		mag[2] = sqrt(mag[0] * mag[1]);

		c = V3DOT(crs[0], crs[1]) / mag[2];
		b11 = -1. * c / mag[0];
		b22 = -1. * c / mag[1];
		b12 = 1. / mag[2];
		V3COPY(e, dr[(k + 1) % NDIM]);
		//
		V3CROSS(crs[2], crs[0], e);
		V3CROSS(crs[3], crs[1], e);
		VSS3ADD(fi[0], crs[2], crs[3], b11, b12);
		VSS3ADD(fi[3], crs[2], crs[3], -1. * b12, -1. * b22);
		//
		V3CROSS(crs[2], dr[3], crs[0]);
		V3CROSS(crs[3], crs[0], dr[4]);
		VSS3ADD(fi[1], crs[2], crs[3], b11, b12);
		V3CROSS(crs[2], dr[3], crs[1]);
		V3CROSS(crs[3], crs[1], dr[4]);
		VSS3ADD(tempDbl, crs[2], crs[3], b12, b22);
		VV3ADD(fi[1], tempDbl);
		V3CROSS(crs[2], crs[0], dr[k]);
		V3CROSS(crs[3], dr[5], crs[0]);
		VSS3ADD(fi[2], crs[2], crs[3], b11, b12);
		V3CROSS(crs[2], crs[1], dr[k]);
		V3CROSS(crs[3], dr[5], crs[1]);
		VSS3ADD(tempDbl, crs[2], crs[3], b12, b22);
		VV3ADD(fi[2], tempDbl);
		for(l = 0; l < 4; l++) {
			VVS3ADD(&P2(memb.f,iMb[ind[l]],0), fi[l], memb.bend.stf);
            if (iMb[ind[l]] < nMbMe) {
                recMb.bendF[iMb[ind[l]]] += V3LEN(fi[l]) * memb.bend.stf;
            }
		}

		V2COPY(&P2A(list.l,list.c,0,2), listInd);
		(list.c)++;
	}
  }
  free(list.l);
}

void CalcMembraneAreaForces(void) {
  int n, k, locInd, *pArr, kind;
  double stf, fac, dr[4][NDIM], xi[NDIM], alpha, lenXi;

  for(n = 0; n < memb.unit.c; n++) {
	pArr = &P2A(memb.unit.l,n,0,3);
	locInd = iMb[pArr[0]];
	kind = ISNUC(memb.idx[locInd]) ? 1 : 0;
	CONT((mbAre.gTgl == 0 && kind == 0) || (nucAre.gTgl == 0 && kind == 1));

	stf = (kind == 1) ? memb.nucArea.stf : memb.area.stf;
	for(k = 0; k < 3; k++) {
		CalcVec(dr[k], &P2(memb.r,iMb[pArr[(k + 1) % NDIM]],0), 
				&P2(memb.r,iMb[pArr[k]],0));
	}
	V3CROSS(xi, dr[2], dr[0]);
	lenXi = V3LEN(xi);
	if ((kind == 0 && mbAre.gTgl == 1) || (kind == 1 && nucAre.gTgl == 1)) { 
		alpha = (0.5 * lenXi - memb.unitEqAr[n]) / memb.unitEqAr[n];
	}
	else {
		alpha = memb.area.val[memb.idx[locInd]];
	}
	fac = -1. * stf * alpha / (2. * lenXi);

	for(k = 0; k < 3; k++) {
		V3CROSS(dr[3], xi, dr[(k + 1) % NDIM]);
		VVS3ADD(&P2(memb.f,iMb[pArr[k]],0), dr[3], fac);
	}
  }
}

/*----------------------------- Calculate force ------------------------------*/

/*---------------------------- Update information ----------------------------*/

void UpdateMembraneBindingSubroutine(int actInd, int mbInd, int side, 
		int sideAct) {
  int locActInd, locMbInd;
  double len;

  locActInd = iAct[actInd];
  locMbInd = iMb[mbInd];
  if (locMbInd > -1) {
    if (side < 0) {
        side = FindElementArray(&P2A(memb.ch,locMbInd,nChMb - nMbAct,nChMb),
                nMbAct, -1, 0, 1);
    }
	if (mbMat.gTgl != 0) {
		P2A(memb.nFA,locMbInd,side,nMbAct) = 1;
	}
	P2A(memb.ch,locMbInd,side + nChMb - nMbAct,nChMb) = actInd;
  }
  if (locActInd > -1) { 
	P2A(act.mbReb,locActInd,sideAct,nRebMb) = mbInd; 
  }
  if (locMbInd > -1 && locActInd > -1) {
	len = CalcDist(&P2(act.r,locActInd,0), &P2(memb.r,locMbInd,0), 0);
	P2A(act.mbRebEq,locActInd,sideAct,nRebMb) = len;
  }

  if (locMbInd > -1 && locMbInd < nMbMe) { (mbReb.cntMe)++; };
  if (((locActInd > -1 && locActInd < nActMe) 
		|| (locMbInd > -1 && locMbInd < nMbMe)) && mpiMethod == 0) { 
	InsertLongChain(mbInd + nAct + nAbp, actInd, minDimDomC * 0.9);
  }
  V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), mbInd, actInd, currTimeStep);
  (noMbDyn.c)++;
}

// Update the binding of actin filaments on membranes.
void UpdateMembraneBinding(void) {
  int m, n, k, side, sideAct, locActInd, locMbInd, CS;
  int *nlPnt, *actCyl, *mbCyl, *pArr;
  double dr[NDIM], len, pReb;
  ListInt chkPair;
  Force *pS;

  nlPnt = neighMb.l;
  chkPair.l = allIntL;
  chkPair.c = 0;
  for(n = 0; n < neighMb.c; n++) {
    if (n > 0) { nlPnt += 2; }
	// Consider the possible bond formation between membrane points and actins
	CONT(!((nlPnt[0] < nAct && nlPnt[1] >= nAct + nAbp) 
			|| (nlPnt[1] < nAct && nlPnt[0] >= nAct + nAbp)));
	side = (nlPnt[0] < nAct) ? 0 : 1;
	CONT(nlPnt[1 - side] >= nAct + nAbp + nUnitMb);
	actCyl = &P2A(act.cyl.l,nlPnt[side],0,2);
	mbCyl = &P2A(memb.unit.l,nlPnt[1 - side] - nAct - nAbp,0,dimMbNuc);
	for(m = 0; m < 2; m++) {
		locActInd = iAct[actCyl[m]];
		CONT(locActInd < 0);
		pArr = &P2A(act.ch,locActInd,0,nChAc);
		// Skip actin monomers
		CONT(pArr[0] < 0 && pArr[1] < 0);
		// Depending on the condition file, the middle part of actin filaments
		// cannot be bound.
		if (mbReb.gTglPa == 0) {
			CONT(pArr[0] > -1 && pArr[1] > -1);
		}
		// Allow the binding only on barbed ends
		CONT(!(pArr[0] < 0 && pArr[1] > -1));
		// If the actin is already bound to membrane, skip it.
		for(k = 0; k < nRebMb; k++) {
			BREAK(P2A(act.mbReb,locActInd,k,nRebMb) < 0);
		}
		CONT(k == nRebMb);
		sideAct = k;
		CS = CheckActinAvailability(actCyl[m], 2);
		CONT(CS > -1);
		for(k = 0; k < dimMbNuc; k++) {
			BREAK(sideAct >= nRebMb);

			locMbInd = iMb[mbCyl[k]];
			// If the information of the membrane point doesn't exist, skip it.
			CONT(locMbInd < 0 || locMbInd >= nMbMe);
			// Check whether the membrane is allowed to have bound actins
			CONT(memb.reb[locMbInd] < 0);
			// Check whether the membrane point has an available binding site
			side = FindElementArray(&P2A(memb.ch,locMbInd,nChMb - nMbAct,nChMb),
					nMbAct, -1, 0, 1);
			CONT(side == -1);
			pS = (ISNUC(memb.idx[locMbInd])) ? &memb.nucSpr2 : &memb.spr2;
			pReb = (ISNUC(memb.idx[locMbInd])) ? nucReb.p : mbReb.p;
			CONT(P2(act.r,locActInd,1) < rGrid[1][0] + dimDom[1] * yLoFA
					|| P2(act.r,locActInd,1) > rGrid[1][0] + dimDom[1] * yHiFA);

			CONT(P2(act.r,locActInd,0) < rGrid[0][0] + dimDom[0] * xLoFA
					|| P2(act.r,locActInd,0) > rGrid[0][0] + dimDom[0] * xHiFA);

			CONT(P2(memb.r,locMbInd,1) < P2(act.r,locActInd,1)); 

			CalcVec(dr, &P2(act.r,locActInd,0), &P2(memb.r,locMbInd,0));
			if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }
			len = V3LEN(dr);
			// If length between actin and membrane point doesn't lie within
			// a certain range, binding cannot occur.
			CONT(len < pS->lo || len >= pS->hi);
			// Check whether this paper was considered before.
			CS = Find2ElementArray(chkPair.l, chkPair.c, actCyl[m], 
					mbCyl[k], 0, 2);
			CONT(CS > -1);
			V2SET(&P2A(chkPair.l,chkPair.c,0,2), actCyl[m], mbCyl[k]);
			(chkPair.c)++;
			// Check a probability
        	CONT(!(genrand_real3() < pReb));
			CS = FindElementArray(noMbDyn.l, noMbDyn.c, mbCyl[k], 0, 3);
			CONT(CS > -1);
			UpdateMembraneBindingSubroutine(actCyl[m], mbCyl[k], side, 
					sideAct);
			V3SET(&P2A(sendMbDyn.l,sendMbDyn.c,0,3), mbCyl[k], actCyl[m], 
					sideAct + 1);
			(sendMbDyn.c)++;

			break;
		}
	}
  }
}

int UpdateMembraneUnbindMatureSubroutine(int actInd, int mbInd, int sideMb,
		int sideAct) {
  int locActInd, locMbInd;
  int k;

  locActInd = iAct[actInd];
  locMbInd = iMb[mbInd];
  if (locMbInd > -1) {
	if (sideMb < 0) {
		sideMb = FindElementArray(&P2A(memb.ch,locMbInd,nChMb - nMbAct,nChMb),
				nMbAct, actInd, 0, 1);
	}
	if (mbMat.gTgl != 0) {
		P2A(memb.nFA,locMbInd,sideMb,nMbAct) = 0;
	}
	P2A(memb.ch,locMbInd,sideMb + nChMb - nMbAct,nChMb) = -1;
  }
  if (locActInd > -1) { 
	if (sideAct < 0) {
		sideAct = FindElementArray(&P2A(act.mbReb,locActInd,0,nRebMb), 
				nRebMb, mbInd, 0, 1);
	}
	for(k = sideAct; k < nRebMb - 1; k++) {
		P2A(act.mbReb,locActInd,k,nRebMb) 
				= P2A(act.mbReb,locActInd,k + 1,nRebMb); 
		P2A(act.mbRebEq,locActInd,k,nRebMb) 
				= P2A(act.mbRebEq,locActInd,k + 1,nRebMb); 
	}
	P2A(act.mbReb,locActInd,nRebMb - 1,nRebMb) = -1;
	P2A(act.mbRebEq,locActInd,nRebMb - 1,nRebMb) = 0.;
  }
  // The pair should be deleted in longCh.l
  if (((locActInd > -1 && locActInd < nActMe) 
			|| (locMbInd > -1 && locMbInd < nMbMe)) && mpiMethod == 0) { 
	DeleteLongChain(mbInd + nAct + nAbp, actInd);
  }
  if (locMbInd > -1 && locMbInd < nMbMe) { (mbUnb.cntMe)++; }
  V3SET(&P2A(noMbDyn.l,noMbDyn.c,0,3), mbInd, actInd, currTimeStep);
  (noMbDyn.c)++;
  return sideAct;
}

// Update the unbinding of actin filaments from membranes.
void UpdateMembraneUnbindMature(void) {
  int n, k, fMag, CS, actInd, locActInd, sideAct;
  double dr[NDIM], pUnb, len, f;
  FOR_MBME(n) {
	CS = FindElementArray(noMbDyn.l, noMbDyn.c, memb.id[n], 0, 3);
	CONT(CS > -1);
	for(k = 0; k < nMbAct; k++) {
		actInd = P2A(memb.ch,n,k + nChMb - nMbAct,nChMb);
		CONT(actInd < 0);
		CS = CheckActinAvailability(actInd, 2);
		CONT(CS > -1);
		locActInd = iAct[actInd];
		// Calculate the unit vector along the chain between membrane point 
	 	// and actin.
	  	CalcVec(dr, &P2(memb.r,n,0), &P2(act.r,locActInd,0));
	 	if (dimMbNuc == 2) { dr[dirNormMbNuc] = 0.; }

		sideAct = FindElementArray(&P2A(act.mbReb,locActInd,0,nRebMb), 
				nRebMb, memb.id[n], 0, 1);
		len = V3LEN(dr);
		f = memb.spr2.stf * (len - P2A(act.mbRebEq,locActInd,sideAct,nRebMb));

		if (f < 0.) { f = 0.; }

		// Maturation
		if (mbMat.gTgl != 0) {
			CS = 1;
			// If actins exist inside the membrane, the maturation can occur
			// only at membrane points fixed in space.
			if (sideMb == 0) {
				if (mbFix.gTgl != 0 || mbDef.gTgl != 0) {
					if (memb.fix[n] < 0) {
						CS = 0;
					}
				}
			}
			if (mbMat.gTglPa == 0 && CS == 1) {
				if (P2A(act.ch,locActInd,0,nChAc) > -1 
						&& P2A(act.ch,locActInd,1,nChAc) > -1) {
					CS = 0;
				}
			}
			if (CS == 1) {
				if (P2A(memb.nFA,n,k,nMbAct) < mbMat.maxNFA) {
					fMag = TrimIntVal((int)f, 0, mbMat.maxF - 1);
					if (genrand_real3() < mbMat.p[fMag]) {
						P2A(memb.nFA,n,k,nMbAct)++;
						mbMat.cntMe2[0]++;
						continue;
					}
				}
			}
		}
		// Unbinding
		if (mbUnb.gTgl != 0) {
			fMag = TrimIntVal((int)f, 0, mbUnb.maxF - 1);
			if (mbMat.gTgl != 0) {
				fMag = (int)(fMag / (P2A(memb.nFA,n,k,nMbAct) + 1));
			}
		    pUnb = mbUnb.p[fMag];
			// Check a probability.
			CONT(!(genrand_real3() < pUnb));
			if (mbMat.gTgl != 0) {
		  		P2A(memb.nFA,n,k,nMbAct)--;
				mbMat.cntMe2[1]++;
				CONT(P2A(memb.nFA,n,k,nMbAct) > 0);
			}
			UpdateMembraneUnbindMatureSubroutine(actInd, memb.id[n], k, 
					sideAct);
			V3SET(&P2A(sendMbDyn.l,sendMbDyn.c,0,3), memb.id[n], actInd, 
					-1 * (sideAct + 1));
			(sendMbDyn.c)++;
		}
	}
  } 
}

void UpdateMembraneUnitArea(void) {
  int n, k, nCh, *pArr;
  double area;

  FOR_MBME(n) {
	pArr = &P2A(memb.ch,n,0,nChMb);
	if (dimMbNuc == 2) { nCh = 2; }
	else { nCh = (pArr[5] < 0) ? 5 : 6; }
	memb.areaL[n] = 0.;
	for(k = 0; k < nCh; k++) {
		if (dimMbNuc == 2) { 
			area = CalcDist(&P2(memb.r,n,0), &P2(memb.r,iMb[pArr[k]],0), 0);
			area *= dimDom[dirNormMbNuc];
			memb.areaL[n] += area;
		}
		else {
			area = CalcTriaArea(&P2(memb.r,n,0), &P2(memb.r,iMb[pArr[k]],0), 
					&P2(memb.r,iMb[pArr[(k + 1) % nCh]],0));
			memb.areaL[n] += area;
		}
	}
	memb.areaL[n] /= (dimMbNuc == 2) ? 2. : 3.;
  }
}

void UpdateMembraneTotalArea(void) {
  int n, k, *pArr;
  double area, *sum, *sumAll, eq;
 
  MALLOC(sum,double,nObjMbNuc);
  MALLOC(sumAll,double,nCpu*nObjMbNuc);
  SetAllValue1dArrayDouble(sum, nObjMbNuc, 0.);

  for(n = 0; n < nMbMe; n++) {
	sum[memb.idx[n]] += memb.areaL[n];
  }
  MPI_Gather(sum, nObjMbNuc, MPI_DOUBLE, sumAll, nObjMbNuc, MPI_DOUBLE, 0, 
		MPI_COMM_WORLD);
  if (rank == 0) {
	for(n = 1; n < nCpu; n++) {
		for(k = 0; k < nObjMbNuc; k++) {
		    sum[k] += P2A(sumAll,n,k,nObjMbNuc);
		}
	}
  }
  MPI_Bcast(sum, nObjMbNuc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (currTimeStep > 1) { 
  	for(n = 0; n < nObjMbNuc; n++) {
		eq = ISNUC(n) ? memb.nucArea.eq : memb.area.eq;
		memb.area.val[n] = (sum[n] - eq) / eq;
	}
  }
  else { 
	memb.area.eq = sum[0]; 
	if (gTglNuc != 0) {
		memb.nucArea.eq = sum[nObjMbNuc / 2]; 
	}
	SetAllValue1dArrayDouble(memb.area.val, nObjMbNuc, 0.);
  }
  free(sum);
  free(sumAll);
}

// Eliminate expired elements in noMbDyn.l
void UpdateNoMembraneDynamicsList(void) {
  int n;

  CheckArraySize(&noMbDyn, &noMbDyn.siz, 3, 1);
  for (n = noMbDyn.c - 1; n >= 0; n--) {
    CONT(!(currTimeStep - P2A(noMbDyn.l,n,2,3) > durNoMbDyn));
	DeleteElementArrayByIndex(noMbDyn.l,&noMbDyn.c,n,3);
  }
}

/*---------------------------- Update information ----------------------------*/

/*------------------------------ Initialization ------------------------------*/


void InitMembraneInformation(void) {
  int m, n, k, cnt, ind, ind2, CS, CS2, nObj, nPerObj, *pArr;
  double r[NDIM], dr[NDIM], v[NDIM], rBnd[NDIM*2], cen[NDIM];
  double ang, por, rad, cNeiEdge;
  // Initialize the position and chain information of membrane points.
  // For 2D membrane
  if (dimMbNuc == 2) {
	if (gTglLoadMbNucData == 0) { 
		nObj = nObjMbNuc / ((gTglNuc != 0) ? 2 : 1);
		for(m = 0; m < nObjMbNuc; m++) {
			ind = m % nObj;
			// Nucleus
			if (ISNUC(m)) {
				ind2 = nObjMbNuc / 2 * nMbPerObj 
						+ (m - nObjMbNuc / 2) * nNucPerObj;
				nPerObj = nNucPerObj;
				rad = radNuc;
			}
			// Membrane
			else {
				ind2 = m * nMbPerObj;
				nPerObj = nMbPerObj;
				rad = radMb;
			}
			ang = 2. * PI / (double)nPerObj;

			if (nObj == 1) {
				V3COPY(&P2(cenMb,m,0), dimDomH);
			}
			else if (nObj == 2) {
				P2(cenMb,m,dirMbNuc[0]) = dimDom[dirMbNuc[0]] 
						* ((ind == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirMbNuc[1]) = dimDom[dirMbNuc[1]] 
						* ((ind == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirNormMbNuc) = dimDomH[dirNormMbNuc];
			}
			else if (nObj == 4) {
				P2(cenMb,m,dirMbNuc[0]) = dimDom[dirMbNuc[0]] 
						* ((ind / 2 == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirMbNuc[1]) = dimDom[dirMbNuc[1]] 
						* ((ind % 2 == 0) ? 0.25 : 0.75);
				P2(cenMb,m,dirNormMbNuc) = dimDomH[dirNormMbNuc];
			}
			for(n = 0; n < nPerObj; n++) {
				P2(rMb,ind2 + n,dirMbNuc[0]) = P2(cenMb,m,dirMbNuc[0])
						+ rad * cos(ang * (double)n);
				P2(rMb,ind2 + n,dirMbNuc[1]) = P2(cenMb,m,dirMbNuc[1])
						+ rad * sin(ang * (double)n);
				P2(rMb,ind2 + n,dirNormMbNuc) = P2(cenMb,m,dirNormMbNuc);
				V2SET(&P2A(chMb,ind2 + n,0,nChMb), ind2 + n + 1, ind2 + n - 1);
				memset(&P2A(chMb,ind2 + n,2,nChMb), -1, sizeof(int) * nMbAct);
				CalcUnitVec(&P2(nDirMb,ind2 + n,0), &P2(cenMb,m,0), 
						&P2(rMb,ind2 + n,0));
				idxMb[ind2 + n] = m;
			}
			P2A(chMb,ind2,1,nChMb) = ind2 + nPerObj - 1;
			P2A(chMb,ind2 + nPerObj - 1,0,nChMb) = ind2;
			for(n = 0; n < nPerObj; n++) {
				V2SET(&P2A(unitMb,ind2 + n,0,2), ind2 + n, 
						ind2 + (n + 1) % nPerObj);
				V3AVG(&P2(nDirUnitMb,ind2 + n,0), &P2(nDirMb,ind2 + n,0), 
						&P2(nDirMb,ind2 + (n + 1) % nPerObj,0));
			}
		}
		nUnitMb = nMb;
	}
  }
  // For 3D membrane
  else {
	InitFlatMembrane();
  }
  // Distribute the information to subdomains.
  memset(iMb, -1, sizeof(int) * nMb);
  cNeiEdge = (motWalk.gTgl != 0) ? neiEdge + 1. : neiEdge;
  V6COPY(rBnd,bnd.r);
  FOR_NDIM(k) {
	CONT(!(pbc[k] == 1));
	if (iCell[k] == 0 && nCell[k] > 1)
	{ P2A(rBnd,0,k,NDIM) = rGrid[k][nGrid[k] - 1]; }
	else if (iCell[k] == nCell[k] - 1 && nCell[k] > 1)
	{ P2A(rBnd,1,k,NDIM) = rGrid[k][0]; }
  }
  nMbMe = 0;
  nMbCp = 0;
  for(m = 0; m < 2; m++) {
	FOR_MB(n) {
		CS = 0;
		CS2 = 0;
		V3COPY(r, &P2(rMb,n,0));
		ExtractConfigSubroutine(r, rBnd, cNeiEdge, m, &CS, &CS2);
		if ((m == 0 && CS == NDIM) ||
				(m == 1 && CS + CS2 == NDIM && CS != NDIM)) {
			ind = (m == 0) ? nMbMe : nMbMe + nMbCp;
			V3COPY(&P2(memb.r,ind,0), &P2(rMb,n,0));
			V3COPY(&P2(memb.nDir,ind,0), &P2(nDirMb,n,0));
			memb.idx[ind] = idxMb[n];
			memb.id[ind] = n;
			iMb[n] = ind;
			if (mbReb.gTgl != 0 || mbUnb.gTgl != 0 
					|| nucReb.gTgl != 0 || nucUnb.gTgl != 0) {
				memb.reb[ind] = rebMb[n]; 
			}
			if (m == 0) { nMbMe++; }
			else { nMbCp++; }
		}
	}
  } 

  FOR_MB(n) {
    CONT(iMb[n] < 0);
    Copy1dArrayInt(&P2A(memb.ch,iMb[n],0,nChMb), &P2A(chMb,n,0,nChMb), nChMb);
	if (dimMbNuc == 3) {
		Copy1dArrayDouble(&P2A(memb.eqLen,iMb[n],0,nChMb - nMbAct), 
				&P2A(eqLenMb,n,0,nChMb - nMbAct), nChMb - nMbAct);
	}
  }
  memb.unit.c = 0;
  memb.unitCp.c = 0;
  for(n = 0; n < nUnitMb; n++) {
 	pArr = &P2A(unitMb,n,0,dimMbNuc);
	for(k = 0; k < dimMbNuc; k++) {
		BREAK(iMb[pArr[k]] > -1 && iMb[pArr[k]] < nMbMe);
	}
	CONT(k == dimMbNuc);
	Copy1dArrayInt(&P2A(memb.unit.l,memb.unit.c,0,dimMbNuc), pArr, dimMbNuc);
	V3COPY(&P2(memb.unitNDir,memb.unit.c,0), &P2(nDirUnitMb,n,0));
	if (dimMbNuc == 3 && (mbAre.gTgl == 1 || nucAre.gTgl == 1)) { 
		memb.unitEqAr[memb.unit.c] = eqArUnitMb[n];
	}
	(memb.unit.c)++;
  }

  radMbInit = radMb;
  radNucInit = radNuc;

  free(rMb);
  free(nDirMb);
  free(chMb);
  free(unitMb);
  free(nDirUnitMb);
  free(idxMb);
  if (dimMbNuc == 3 && (mbAre.gTgl == 1 || nucAre.gTgl == 1)) { 
	free(eqArUnitMb);
  }
  if (mbReb.gTgl != 0) { free(rebMb); }
  if (dimMbNuc == 3) {
	free(eqLenMb);
  }
  UpdateSubdomSectionLocation(0);
}

void Init3dMembraneInformation(void){
  // Array used for initializing the vertices of the icosahedron
  int initTriaArr[60] = {0,2,1,0,3,2,0,4,3,0,5,4,0,1,5,1,2,7,2,3,8,3,4,9,4,5,
		10,5,1,6,1,7,6,2,8,7,3,9,8,4,10,9,5,6,10,6,7,11,7,8,11,8,9,11,9,10,11,
		10,6,11};
  int mm, m, n, k, l, CS, nMb2, nUnitMb2, *pArr, *pArr2, nCh, nObj, nPerObj;
  int nUnitMbPerObj, nUnitNucPerObj, nUnitPerObj;
  int ind[2], vert[6], order[12] = {0,3,5,3,1,4,5,4,2,5,3,4};
  int *chMb2, *unitMb2;	
  double theta, sinThe, cosThe, phi, len, lev, rad, min;
  double dr[4][NDIM], r[NDIM], rIcos[12][NDIM], ang[2];
  double *rMb2, *nDirMb2, *eqLenMb2, *nDirUnitMb2, *eqArUnitMb2;
  ListInt newVert;

  nUnitMbPerObj = 20 * (int)pow(4., levMb);
  nUnitNucPerObj = 20 * (int)pow(4., levNuc);
  nObj = nObjMbNuc / ((gTglNuc != 0) ? 2 : 1);
  for(mm = 0; mm < ((gTglNuc != 0) ? 2 : 1); mm++) {
	if (mm == 0) {
		lev = levMb;
		nPerObj = nMbPerObj;
		rad = radMb;
		nUnitPerObj = nUnitMbPerObj;
	}
	else {
		lev = levNuc;
		nPerObj = nNucPerObj;
		rad = radNuc;
		nUnitPerObj = nUnitNucPerObj;
	}
	MALLOC(rMb2,double,NDIM*nPerObj);
	MALLOC(nDirMb2,double,NDIM*nPerObj);
	MALLOC(chMb2,int,6*nPerObj);
	MALLOC(eqLenMb2,double,6*nPerObj);
	MALLOC(unitMb2,int,3*nUnitPerObj);
	MALLOC(nDirUnitMb2,double,NDIM*nUnitPerObj);
	if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
		MALLOC(eqArUnitMb2,double,nUnitPerObj);
	}

	memset(chMb2, -1, (sizeof(int) * 6 * nPerObj));

	theta = DEG2RAD(26.56505117707799);
	sinThe = sin(theta);
	cosThe = cos(theta);
	// The lower and upper vertex of icosahedron
	V3SET(rIcos[0], 0, 0, -1.); 
	V3SET(rIcos[11], 0, 0, 1.); 	
	// Lower pentagon of icosahedron
	phi = PI / 5;
	for (n = 1; n < 6; n++) {
		V3SET(rIcos[n], cosThe * cos(phi), cosThe * sin(phi), -1 * sinThe);
		phi += 0.4 * PI;
	}
	// Upper pentagon of icosahedron
	phi = 0;
	for (n = 6; n < 11; n++) {
	  	V3SET(rIcos[n], cosThe * cos(phi), cosThe * sin(phi), sinThe);
	  	phi += 0.4 * PI;
	}
	for (n = 0; n < 12; n++) {
		V3SCALE(rIcos[n], rad);
		V3COPY(&P2(rMb2,n,0), rIcos[n]);
	}
	Copy1dArrayInt(unitMb2, initTriaArr, 60);

	nMb2 = 12;
	nUnitMb2 = 20;

	for(m = 0; m < lev; m++) {
		MALLOC(newVert.l,int,(int)(30 * pow(4, m + 1)) * 3);
		newVert.c = 0;
		for(n = 0; n < nUnitMb2; n++) { 
			V3COPY(vert, &P2A(unitMb2,n,0,3));
			// Determines the next 3 nodes for the the current triangle's 
			// subdivision
			for(k = 0; k < 3; k++){
				V2SET(ind, k, (k + 1) % 3);
				for(l = 0; l < newVert.c; l++) {
					pArr = &P2A(newVert.l,l,0,3);
					BREAK((pArr[0] == vert[ind[0]] && pArr[1] == vert[ind[1]])
							|| (pArr[0] == vert[ind[1]] 
							&& pArr[1] == vert[ind[0]]));
				}
				if (l < newVert.c) {
					vert[k + 3] = P2A(newVert.l,l,2,3);
				}
				else {
					vert[k + 3] = nMb2;
					V3ADD(r, &P2(rMb2,vert[ind[0]],0), 
							&P2(rMb2,vert[ind[1]],0));
					NormVec(r);
					V3SCALE(r, rad);
					V3COPY(&P2(rMb2,nMb2,0), r);
					V3SET(&P2A(newVert.l,newVert.c,0,3), vert[ind[0]], 
							vert[ind[1]], nMb2);
					(newVert.c)++;
					nMb2++;
				}
			}
			for(k = nUnitMb2 - 1; k >= n; k--){
				V3COPY(&P2A(unitMb2,k + 3,0,3), &P2A(unitMb2,k,0,3));
			}		
			for(k = 0; k < 12; k++){
				unitMb2[n * 3 + k] = vert[order[k]];
			}		
			nUnitMb2 += 3;
			n += 3;
		}
		free(newVert.l);
	}
	// Error check
	if (nPerObj != nMb2 || nUnitPerObj != nUnitMb2) {
		printf("Error: %d vs %d\n", nMb, nPerObj);
	}

	// Offset the membrane position to the center of a domain
	for(n = 0; n < nMb2; n++) {
		VS3COPY(&P2(nDirMb2,n,0), &P2(rMb2,n,0), -1.);
		NormVec(&P2(nDirMb2,n,0));
	}

	for(m = 0; m < nUnitMb2; m++) {
		qsort(&P2A(unitMb2,m,0,3), 3, sizeof(int), CompInt);
	    pArr = &P2A(unitMb2,m,0,3);
		CalcVec(dr[0], &P2(rMb2,pArr[1],0), &P2(rMb2,pArr[0],0));
		CalcVec(dr[1], &P2(rMb2,pArr[2],0), &P2(rMb2,pArr[0],0));
		V3CROSS(dr[2], dr[0], dr[1]);	
		if (V3COS(dr[2], &P2(nDirMb2,pArr[0],0)) < 0.) {
			SwapInt(&pArr[1], &pArr[2]);
			V3REVSIGN(dr[2]);
		}
		if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
			eqArUnitMb2[m] = 0.5 * V3LEN(dr[2]);
		}
		NormVec(dr[2]);
		V3COPY(&P2(nDirUnitMb2,m,0), dr[2]);
		
	    for(n = 0; n < 3; n++) {
	        ind[0] = pArr[n];
	        ind[1] = pArr[(n + 1) % 3];
	        len = CalcDist(&P2(rMb2,ind[0],0), &P2(rMb2,ind[1],0), 0);
	        for(k = 0; k < 2; k++) {
	            pArr2 = &P2A(chMb2,ind[k],0,6);
	            CS = -1;
	            for(l = 0; l < 6; l++) {
	                if (pArr2[l] == ind[1 - k]) {
	                    CS = 1;
	                    break;
	                }
	                BREAK(pArr2[l] == -1);
	            }
	            if (CS == -1 && l < 6) {
	                pArr2[l] = ind[1 - k];
	                P2A(eqLenMb2,ind[k],l,6) = len;
	            }
	        }
	    }
	}
 
	for(m = 0; m < nMb2; m++) {
		nCh = (P2A(chMb2,m,5,6) < 0) ? 5 : 6;
		CalcVec(dr[3], &P2(rMb2,P2A(chMb2,m,0,6),0), &P2(rMb2,m,0));
		V3CROSS(dr[0], dr[3], &P2(nDirMb2,m,0));
		for(n = 1; n < nCh; n++) {
			for(k = n + 1; k < nCh; k++) {
				for(l = 0; l < 2; l++) {
					CalcVec(dr[3], &P2(rMb2,P2A(chMb2,m,(l == 0) ? n : k,6),0), 
							&P2(rMb2,m,0));
					V3CROSS(dr[l + 1], dr[3], &P2(nDirMb2,m,0));
					ang[l] = V3ANG(dr[0], dr[l + 1]);
					V3CROSS(dr[3], dr[0], dr[l + 1]);
					if (V3DOT(&P2(nDirMb2,m,0), dr[3]) > 0.) {
						ang[l] = 2 * PI - ang[l];
					}
				}
				if (ang[0] > ang[1]) {
				 	SwapInt(&P2A(chMb2,m,n,6), &P2A(chMb2,m,k,6));
				 	SwapDbl(&P2A(eqLenMb2,m,n,6), &P2A(eqLenMb2,m,k,6));
				}	
			}
		}
	}
	for(m = 0; m < nObj; m++) {
		if (nObj == 1) {
			V3COPY(cenMb, dimDomH);
		}
		else if (nObj == 2) {
			VS3COPY(&P2(cenMb,m,0), dimDom, (m % nObj == 0) ? 0.25 : 0.75);
		}
		ind[0] = mm * nMbPerObj * nObjMbNuc / 2 + m * nPerObj;
		for(n = 0; n < nPerObj; n++) {
			V3ADD(&P2(rMb,ind[0] + n,0), &P2(rMb2,n,0), &P2(cenMb,m,0));
			V3COPY(&P2(nDirMb,ind[0] + n,0), &P2(nDirMb2,n,0));
			for(k = 0; k < 6; k++) {
				P2A(eqLenMb,ind[0] + n,k,6) = P2A(eqLenMb2,n,k,6);
				P2A(chMb,ind[0] + n,k,nChMb) = P2A(chMb2,n,k,6) 
						+ (P2A(chMb2,n,k,6) > -1 ? ind[0] : 0);
			}
	  		memset(&P2A(chMb,ind[0] + n,6,nChMb), -1, 
					sizeof(int) * (nChMb - 6));
			idxMb[ind[0] + n] = m + mm * nObj;
		}
		ind[1] = mm * nUnitMbPerObj * nObjMbNuc / 2 + m * nUnitPerObj;
		for(n = 0; n < nUnitPerObj; n++) {
			for(k = 0; k < 3; k++) {
				P2A(unitMb,ind[1] + n,k,3) = P2A(unitMb2,n,k,3) + ind[0];
			}
			V3COPY(&P2(nDirUnitMb,ind[1] + n,0), &P2(nDirUnitMb2,n,0));
			if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
				eqArUnitMb[ind[1] + n] = eqArUnitMb2[n];
			}
		}
	}

	min = POS_LARGE_VALUE;
	for(n = 0; n < nMb; n++) {
		for(k = 0; k < 6; k++) {
			CONT(P2A(chMb,n,k,nChMb) < 0);
			CONT(P2A(eqLenMb,n,k,6) > min);
			min = P2A(eqLenMb,n,k,6);
		}
	}
	if (0.5 * min < memb.thk) {
		if (rank == 0) {
			printf("Error: the minimum length of membrane chain (%g) should be "
					"sufficiently larger than its thickness (%g)!\n", 
					min, memb.thk);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(0);
	}

	free(rMb2);
	free(nDirMb2);
	free(eqLenMb2);
	free(chMb2);
	free(unitMb2);
	free(nDirUnitMb2);
	if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
		free(eqArUnitMb2);
	}
  }
}

/*------------------------------ Initialization ------------------------------*/

/*----------------------- Process transfered messages ------------------------*/

void UpdateMembraneDynamicsEventsSubroutine(ListInt *all) {
  int n, side, CS, CS2, mode;
  int mbInd, actInd, locMbInd, locActInd, mbInd2;
  int actRankMol, mbRankMol, mbRankMol2; 
  ListInt confMb;

  confMb.l = allIntL;
  confMb.c = 0;
  for(n = sendMbDyn.c; n < all->c; n++) {
	// mode = 0: unbind, 2: bind
	mode = (P2A(all->l,n,2,3) < 0) ? 0 : 1;
	mbInd = P2A(all->l,n,0,3);
	actInd = P2A(all->l,n,1,3);
	side = abs(P2A(all->l,n,2,3)) - 1;	
	locMbInd = iMb[mbInd];
	locActInd = iAct[actInd];
	if (locMbInd >= nMbMe || locMbInd < 0) {
		CS = 1;
		if (locActInd > -1) {
			actRankMol = CalcRankMolecule(&P2(act.r,locActInd,0));
			mbInd2 = P2A(act.mbReb,locActInd,side,nRebMb);
			if (mbInd2 > -1 && mbInd2 != mbInd) {
				CS2 = Find2ElementArray(confMb.l, confMb.c, locActInd, 
						side, 0, 2);
				if (CS2 == -1) {
					V2SET(&P2A(confMb.l,confMb.c,0,2), locActInd, side);
					(confMb.c)++;
				}
				if (iMb[mbInd2] > -1) {
					mbRankMol2 = CalcRankMolecule(&P2(memb.r,iMb[mbInd2],0));
					if (mbRankMol2 != actRankMol) {
						UpdateMembraneUnbindMatureSubroutine(actInd, 
								mbInd2, -1, side);
						if (locMbInd > -1) {
							mbRankMol 
									= CalcRankMolecule(&P2(memb.r,locMbInd,0));
							CS = (actRankMol == mbRankMol) ? 1 : 0;
						}
						else { CS = 1; }
					}
					else { CS = 0; }
				}
			}
			else if (mbInd2 < 0) {
				CS2 = Find2ElementArray(confMb.l, confMb.c, locActInd, 
						side, 0, 2);
				if (CS2 > -1) {
					if (locMbInd > -1) {
						mbRankMol = CalcRankMolecule(&P2(memb.r,locMbInd,0));
						CS = (actRankMol == mbRankMol) ? 1 : 0;
					}
					else { CS = 1; }
				}
				else { CS = 1; }
			}
			// Modify act.ch
			if (CS == 1) {
				if (mode == 0) {
					UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, 
							-1, side);
				}
				else {
					UpdateMembraneBindingSubroutine(actInd, mbInd, -1, side);
				}
				if (locMbInd >= nMbMe) { CS = 2; }
			}
		}
		if (locMbInd >= nMbMe && CS == 1) {
			if (mode == 0) {
				UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, -1, side);
			}
			else {
				UpdateMembraneBindingSubroutine(actInd, mbInd, -1, side);
			}
		}
	}
	if (locMbInd > -1 && locMbInd < nMbMe) {
		if (mode == 0) {
			UpdateMembraneUnbindMatureSubroutine(actInd, mbInd, -1, side);
		}
	}
  }
  sendMbDyn.c = 0;
}
void UpdateMembraneDynamicsEvents(void) {
  int n, sizeArr;
  ListInt all;

  sizeArr = sendMbDyn.siz;
  sizeArr *= (mpiMethod == 0) ? cntAdjRank[0] : 2;
  MALLOC(all.l,int,sizeArr);

  all.c = sendMbDyn.c;
  for(n = 0; n < sendMbDyn.c; n++) {
	V3COPY(&P2A(all.l,n,0,3), &P2A(sendMbDyn.l,n,0,3));
  }
  CollectArrayIntFromAdjacentSubdomain(&all, 3);
  UpdateMembraneDynamicsEventsSubroutine(&all);

  free(all.l);
}

/*----------------------- Process transfered messages ------------------------*/

/*----------------- Interactions between membrane and others -----------------*/
double HowManyInOutMembrane(double *rBnd) {
  int k, cnt, nObj;
  int ind[NDIM], begin[NDIM], end[NDIM];
  double dist, thres[2], r[NDIM], resol, min, max, bound, cenDom[NDIM];

  resol = 0.5;
  if (bnd.gTglRnd != 0) {
	FOR_NDIM(k) {
    	cenDom[k] = 0.5 * (rGrid[k][0] + rGrid[k][nGrid[k] - 1]);
	}
  }
  if (sideMb == 0) {
	thres[0] = (radMb - 0.5 * memb.thk) * (1. - thkMbActNuc);
	thres[1] = radMb - 0.5 * memb.thk;
	if (gTglNuc != 0) {
		if (thres[0] < radNuc + 0.5 * memb.nucThk) {
			thres[0] = radNuc + 0.5 * memb.nucThk;
		}
	}
  }
  else {
	thres[0] = radMb + 0.5 * memb.thk;
	thres[1] = POS_LARGE_VALUE;
  }
  nObj = nObjMbNuc / (gTglNuc != 0 ? 2 : 1);

  if (sideMb == 0) {
  min = POS_LARGE_VALUE;
  max = NEG_LARGE_VALUE;
  for(k = 0; k < nObj; k++) {
	min = (min > P2(cenMb,k,0) - radMb) ? P2(cenMb,k,0) - radMb : min;
	max = (max < P2(cenMb,k,0) + radMb) ? P2(cenMb,k,0) + radMb : max;
  }

  for(k = 0; k < dimMbNuc; k++) {
	bound = (min > P2A(rBnd,0,dirMbNuc[k],NDIM)) 
			? min : P2A(rBnd,0,dirMbNuc[k],NDIM);
	begin[k] = (int)ceil(bound / resol);
	bound = (max < P2A(rBnd,1,dirMbNuc[k],NDIM)) 
			? max : P2A(rBnd,1,dirMbNuc[k],NDIM);
	end[k] = (int)floor(bound / resol);
  }
  }
  else {
	for(k = 0; k < dimMbNuc; k++) {
		begin[k] = (int)ceil(P2A(rBnd,0,dirMbNuc[k],NDIM) / resol);
		end[k] = (int)floor(P2A(rBnd,1,dirMbNuc[k],NDIM) / resol);
	}
  }

  cnt = 0;
  if (dimMbNuc == 3) {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			for(ind[2] = begin[2]; ind[2] < end[2]; ind[2]++) {	
				VS3COPY(r, ind, resol);
				if (bnd.gTglRnd != 0) {
					dist = CalcDist(cenDom, r, 0);
					CONT(dist > bnd.radRnd);
				}
				for(k = 0; k < nObj; k++) {
					dist = CalcDist(&P2(cenMb,k,0), r, 0);
					if (dist >= thres[0] && dist <= thres[1]) { 
						if (sideMb == 0) {
							cnt++; 
							break;
						}
						else { continue; }
					}
					if (sideMb == 1) { 
						cnt++;
						break;
					}
				}
			}
		}
	}
	return (cnt * CUBE(resol));
  }
  else {
	for(ind[0] = begin[0]; ind[0] < end[0]; ind[0]++) {
		for(ind[1] = begin[1]; ind[1] < end[1]; ind[1]++) {
			r[dirMbNuc[0]] = ind[0] * resol;
			r[dirMbNuc[1]] = ind[1] * resol;
			r[dirNormMbNuc] = cenMb[dirNormMbNuc];
			if (bnd.gTglRnd != 0) {
				dist = CalcDist(cenDom, r, 0);
				CONT(dist > bnd.radRnd);
			}
			for(k = 0; k < nObj; k++) {
				dist = CalcDist(&P2(cenMb,k,0), r, 0);
				if (dist >= thres[0] && dist <= thres[1]) { 
					if (sideMb == 0) {
						cnt++; 
						break;
					}
					else { continue; }
				}
				if (sideMb == 1) { 
					break;
				}
			}
			if (sideMb == 1 && k == nObj) {
				cnt++;
			}
		}
	}
	return (cnt * SQR(resol));
  }
}

/*----------------- Interactions between membrane and others -----------------*/

/*-------------------------------- Check errors ------------------------------*/

void CheckMembraneError(int mode) {
  int n, k, ind, CS, *pArr;

  FOR_MBME(n) {
    if ((isnan(P2(memb.r,n,0)) != 0 || isnan(P2(memb.r,n,1)) != 0
            || isnan(P2(memb.r,n,2)) != 0) 
			|| (fabs(P2(memb.f,n,0)) > 2e5 || fabs(P2(memb.f,n,1)) > 2e5 
			|| fabs(P2(memb.f,n,2)) > 2e5) 
			|| isnan(P2(memb.f,n,0)) != 0 || isnan(P2(memb.f,n,1)) != 0  
			|| isnan(P2(memb.f,n,2)) != 0 
    		|| isnan(P2(memb.nDir,n,0)) != 0 || isnan(P2(memb.nDir,n,1)) != 0
            || isnan(P2(memb.nDir,n,2)) != 0) {
//			|| P2A(memb.ch,n,0,nChMb) < 0 || P2A(memb.ch,n,1,nChMb) < 0) {
		printf("currTimeStep = %lld, mode = %d, rank = %d\n", currTimeStep, mode, rank);
		RecordErrorElement(memb.id[n], 2);
	    printf("\n");
    	exit(-1);
	}
  }

  for(n = 0; n < memb.unit.c; n++) {
    pArr = &P2A(memb.unit.l,n,0,dimMbNuc);
	if (iMb[pArr[0]] < 0 || iMb[pArr[0]] >= nMbMe + nMbCp 
			|| iMb[pArr[1]] < 0 || iMb[pArr[1]] >= nMbMe + nMbCp) {
		printf("rank = %d, currTimeStep = %lld\n", rank, currTimeStep);
		printf("UNIT: %d\n", mode);
		for(k = 0; k < 2; k++) {
			RecordErrorElement(pArr[k], 2);
		}
	}
  }

  for(n = 0; n < nActMe; n++) {
	for(k = 0; k < nRebMb; k++) {
		ind = P2A(act.mbReb,n,k,nRebMb);
		CONT(ind < 0);
		CONT(!(iMb[ind] > -1 && iMb[ind] < nMbMe));
		CS = FindElementArray(&P2A(memb.ch,iMb[ind],dimMbNuc,nChMb), 1, 
				act.id[n], 0, 1);
		CONT(!(CS == -1));
		printf("mode = %d\n", mode);
		printf("act: %d (%d/%d)\n", act.id[n], n, nActMe);
		RecordErrorElement(ind, 2);
		exit(-1);
	}
  }
}

void InitFlatMembrane(void) {
  int i, j, k, iPrev, iNext, jNext, iPrevNext, kadj, cnt, nBd, ind;
  int CS, CS2, CS3, CS4;
  double xdist, ydist, dist, dr[3][3];

  nBd = nChMb - nMbAct;
  ydist = dimDom[1] / (double)memb.nXY[1];
  xdist = dimDom[0] / (double)memb.nXY[0];
  memset(chMb, -1, sizeof(int) * nMb * nChMb);

  for (j = 0; j < memb.nXY[1]; j++) {
	for (i = 0; i < memb.nXY[0]; i++) {
		k = memb.nXY[0] * j + i;
		P2(rMb,k,0) = i * xdist + (((j % 2) == 0) ? 0 : 0.5) * xdist;
		P2(rMb,k,1) = j * ydist;
		P2(rMb,k,2) = 0.;
	}
  }
  // Assign chain information between matrix points
  for(j = 0; j < memb.nXY[1]; j++) {
	CS = 1;
	jNext = j + 1;
	if (jNext == memb.nXY[1]) {
		if (pbc[1] == 0) { CS = -1; }
		else { jNext = 0; }
	}
	for(i = 0; i < memb.nXY[0]; i++) {
 		k = memb.nXY[0] * j + i;

		CS2 = 1;
		iNext = i + 1;
		if (iNext == memb.nXY[0]) {
			if (pbc[0] == 0) { CS2 = -1; }
			else { iNext = 0; }
		}

		CS3 = 1;
		CS4 = 1;
		if (j % 2 == 1) { 
			iPrev = i; 
			iPrevNext = i + 1;
			if (iPrevNext == memb.nXY[0]) {
				if (pbc[0] == 0) { CS4 = -1; }
				else { iPrevNext = 0; }
			}
		}
		else {
			iPrev = i - 1;
			if (iPrev == -1) {
				if (pbc[0] == 0) { CS3 = -1; }
				else { iPrev = memb.nXY[0] - 1; }
			}
			iPrevNext = i;
		}
		if (CS == 1 && CS3 == 1) {
			kadj = jNext * memb.nXY[0] + iPrev;
			dist = CalcDist(&P2(rMb,k,0), &P2(rMb,kadj,0), 0);
			P2A(eqLenMb,k,0,nBd) = dist;
			P2A(eqLenMb,kadj,3,nBd) = dist;
			P2A(chMb,k,0,nChMb) = kadj;
			P2A(chMb,kadj,3,nChMb) = k;
		}
		if (CS == 1 && CS4 == 1) {
			kadj = jNext * memb.nXY[0] + iPrevNext;
			dist = CalcDist(&P2(rMb,k,0), &P2(rMb,kadj,0), 0);
			P2A(eqLenMb,k,1,nBd) = dist;
			P2A(eqLenMb,kadj,4,nBd) = dist;
			P2A(chMb,k,1,nChMb) = kadj;
			P2A(chMb,kadj,4,nChMb) = k;
		}
		if (CS2 == 1) {
			kadj = j * memb.nXY[0] + iNext;
			dist = CalcDist(&P2(rMb,k,0), &P2(rMb,kadj,0), 0);
			P2A(eqLenMb,k,2,nBd) = dist;
			P2A(eqLenMb,kadj,5,nBd) = dist;
			P2A(chMb,k,2,nChMb) = kadj;
			P2A(chMb,kadj,5,nChMb) = k;
		}
	}
  }
  for(i = 0; i < nMb; i++) {
	memset(&P2A(chMb,i,6,nChMb), -1, sizeof(int) * nMbAct);
	idxMb[i] = 0;
	V3SET(&P2(nDirMb,i,0), 0., 0., 1.);
  }
  cnt = 0;
  for(i = 0; i < nMb; i++) {
	for(j = 2; j < 4; j++) {
		CONT(P2A(chMb,i,j,nChMb) < 0 || P2A(chMb,i,j + 1,nChMb) < 0);
		V3SET(&P2A(unitMb,cnt,0,3), i, P2A(chMb,i,j,nChMb), 
				P2A(chMb,i,j + 1,nChMb));
		V3SET(&P2(nDirUnitMb,cnt,0), 0., 0., 1.);
		cnt++;
		if (mbAre.gTgl == 1 || nucAre.gTgl == 1) { 
			CalcVec(dr[0], &P2(rMb,i,0), &P2(rMb,P2A(chMb,i,j,nChMb),0));
			CalcVec(dr[1], &P2(rMb,i,0), &P2(rMb,P2A(chMb,i,j + 1,nChMb),0));
			V3CROSS(dr[2], dr[0], dr[1]);	
			eqArUnitMb[cnt] = 0.5 * V3LEN(dr[2]);
		}
	}
  }
  if (cnt != nUnitMb) {
	if (rank == 0) {
		printf("Error: %d %d\n", cnt, nUnitMb);
	}
	exit(-1);
  }
	
  if (mbReb.gTgl != 0 && gTglLoadMbNucData == 0) {
	if (rank == 0) {
		memset(rebMb, -1, sizeof(int) * nMb);
		cnt = 0;
		while(cnt < (int)(mbReb.por * nMb)) {
			ind = GenRandIntIndex(nMb);
			if (rebMb[ind] == -1) {
				rebMb[ind] = 1;
				cnt++;
			}
		}
	}
	MPI_Bcast(rebMb, nMb, MPI_INT, 0, MPI_COMM_WORLD);
  }
  V3SET(cenMb, dimDom[0] * 0.5, dimDom[1] * 0.5, 0);
  
}

void CheckMembraneErrors(int mode) {
  int n, k, actInd, locActInd, sideAct;
  for(n = 0; n < nMbMe; n++) {
	for(k = 0; k < nMbAct; k++) {
		CONT(P2A(memb.ch,n,nChMb - nMbAct+k,nChMb) < 0);
		actInd = P2A(memb.ch,n,nChMb - nMbAct+k,nChMb);
		locActInd = iAct[actInd];
		sideAct = FindElementArray(&P2A(act.mbReb,locActInd,0,nRebMb), 
				nRebMb, memb.id[n], 0, 1);
	}
  }
}
