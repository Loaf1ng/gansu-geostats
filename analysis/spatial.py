# -*- coding: utf-8 -*-
"""Training-only fitting, geographic holdouts, nested selection and sensitivities.
Each fold refits the variogram and preprocessing on its training stations.
"""
import os,json,csv,warnings
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist,squareform,cdist
from scipy.cluster.vq import kmeans2
from itertools import combinations
from scipy.optimize import least_squares, minimize_scalar
import kriging as g
ROOT = g.ROOT
SITES, XY, N = g.SITES, g.XY.copy(), g.N
PREV = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
SPECS=[{'id':r['model'],'variable':r['variable'],'kind':r['variogram'],'transform':r['transform'],'trend':r['trend'],'k':r['neighbors'],'co':r['method']=='CK'} for r in PREV['model_metrics']]
def dump(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def csvout(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def groups(x,k,seed=20261030):
    best=None
    for start in range(20):
        try:
            center,label=kmeans2(x,k,iter=100,minit='++',missing='raise',rng=np.random.default_rng(seed+start))
            loss=float(np.sum((x-center[label])**2))
            if best is None or loss<best[0]:best=(loss,center,label)
        except Exception:continue
    assert best is not None
    center,label=best[1:];order=np.lexsort((center[:,1],center[:,0]));mapping=np.argsort(order)
    return mapping[label]
OPTIONS={'range_multiplier':2.,'bins':10}
def fit_variogram(u,kind,aux=None):
    pairs=g.PAIR_DIST;ii,jj=g.II,g.JJ;du=u[ii]-u[jj];da=None if aux is None else aux[ii]-aux[jj]
    edges=np.linspace(0,pairs.max()/2,OPTIONS['bins']+1);emp=[]
    for a,b in zip(edges[:-1],edges[1:]):
        use=(pairs>=a)&(pairs<=b if b==edges[-1] else pairs<b)
        if use.sum()<5:continue
        emp.append({'h':float(pairs[use].mean()),'count':int(use.sum()),'g1':float((.5*du[use]**2).mean()),'g2':None if da is None else float((.5*da[use]**2).mean()),'g12':None if da is None else float((.5*du[use]*da[use]).mean())})
    if len(emp)<3:raise ValueError('Fewer than three supported lag bins; no fallback/imputation')
    h=np.array([r['h'] for r in emp]);weight=np.sqrt([r['count'] for r in emp]);weight/=weight.max();v1=np.array([r['g1'] for r in emp])
    amin=max(5.,pairs.min()/2);amax=float(pairs.max()*OPTIONS['range_multiplier'])
    if aux is None:
        lo=[amin,1e-6,1e-8];hi=[amax,10.,3.]
        def loss(p):return weight*(p[2]+p[1]*(1-g.corr(h,kind,p[0]))-v1)
    else:
        v2=np.array([r['g2'] for r in emp]);v12=np.array([r['g12'] for r in emp]);lo=[amin,1e-6,1e-6,1e-8,1e-8,-.995];hi=[amax,10.,10.,3.,3.,.995]
        def loss(p):
            sh=1-g.corr(h,kind,p[0]);cross=p[5]*np.sqrt(p[1]*p[2])
            return np.r_[weight*(p[3]+p[1]*sh-v1),weight*(p[4]+p[2]*sh-v2),weight*(cross*sh-v12)]
    fits=[]
    for factor in [.2,.6,1.2]:
        r=float(np.clip(pairs.max()*factor,amin*1.01,amax*.99));initial=[r,1.,.05] if aux is None else [r,1.,1.,.05,.05,float(np.clip(np.corrcoef(u,aux)[0,1],-.8,.8))]
        fitted=least_squares(loss,initial,bounds=(lo,hi),max_nfev=3000,x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10)
        if fitted.success and np.isfinite(fitted.x).all():fits.append(fitted)
    if not fits:raise ValueError('No converged variogram fit')
    best=min(fits,key=lambda r:sum(r.fun**2));p=best.x
    result={'range_km':float(p[0]),'range_upper_km':amax,'range_upper_hit':bool(p[0]>=amax*.999),'partial1':float(p[1]),'nugget1':float(p[2] if aux is None else p[3]),'weighted_sse':float(sum(best.fun**2)),'supported_bins':len(emp),'bound_hits':[i for i in range(len(p)) if np.isclose(p[i],lo[i],rtol=1e-5,atol=1e-7) or np.isclose(p[i],hi[i],rtol=1e-5,atol=1e-7)]}
    if aux is not None:result.update(partial2=float(p[2]),nugget2=float(p[4]),rho=float(p[5]),cross_partial=float(p[5]*np.sqrt(p[1]*p[2])))
    return result,emp
g.fit_variogram=fit_variogram
def configure(train,options=None):
    global OPTIONS
    OPTIONS={'range_multiplier':2.,'bins':10};OPTIONS.update(options or {})
    g.SITES=[SITES[i] for i in train];g.N=len(train);g.XY=XY[train].copy();g.ORIGIN=g.XY.mean(axis=0);g.AXIS_SCALE=g.XY.std(axis=0,ddof=1)
    g.II,g.JJ=np.triu_indices(g.N,1);g.PAIR_DIST=pdist(g.XY);g.D=squareform(g.PAIR_DIST);g.MAXLAG=float(g.PAIR_DIST.max()/2);g.EDGES=np.linspace(0,g.MAXLAG,OPTIONS['bins']+1);g.FIT_CACHE={}
    for v in ['precip','evap']:
        z=np.array([p[v] for p in g.SITES]);g.ESDA['variables'][v]['boxcox_lambda_mle']=float(stats.boxcox(z)[1])
def model(spec):return g.Model(spec['variable'],kind=spec['kind'],transform=spec['transform'],trend=spec['trend'],k=spec['k'],co=spec['co'])
def predict(spec,train,test):
    assert not set(train)&set(test)
    m=model(spec);pred,se=m.predict(XY[test]);obs=np.array([SITES[i][spec['variable']] for i in test]);rows=[]
    for i,o,p,e in zip(test,obs,pred,se):
        valid=bool(np.isfinite(p) and np.isfinite(e) and e>0)
        rows.append({'station_id':int(SITES[i]['id']),'station':SITES[i]['name'],'observed':float(o),'predicted':float(p) if valid else None,'se':float(e) if valid else None,'error':float(p-o) if valid else None,'valid':valid})
    return rows,m
def metrics(rows):
    valid=[r for r in rows if r['valid']];complete=len(valid)==N
    err=np.array([r['error'] for r in valid]);se=np.array([r['se'] for r in valid])
    return {'n_valid':len(valid),'n_expected':N,'complete':complete,'RMSE':float(np.sqrt(np.mean(err**2))) if complete else None,'ME':float(err.mean()) if complete else None,'ASE':float(se.mean()) if complete else None,'RMSSE':float(np.sqrt(np.mean((err/se)**2))) if complete else None,'partial_valid_RMSE':float(np.sqrt(np.mean(err**2))) if len(valid) else None}
def main():
    fold5=groups(XY,5);fold8=groups(XY,8)
    csvout(ROOT/'tables'/'station_blocks.csv',[{'id':int(p['id']),'station':p['name'],'x_km':float(XY[i,0]),'y_km':float(XY[i,1]),'block5':int(fold5[i])+1,'block8':int(fold8[i])+1} for i,p in enumerate(SITES)])
    schemes={'refit_LOO':np.arange(N),'spatial5':fold5,'spatial8':fold8,'spatial5_buffer100km':fold5}
    allrows=[];allmetrics=[];foldrecords=[];fitrecords=[]
    for scheme,labels in schemes.items():
        scheme_rows={s['id']:[] for s in SPECS}
        for fold in np.unique(labels):
            test=np.flatnonzero(labels==fold);train=np.flatnonzero(labels!=fold)
            if 'buffer' in scheme:train=train[cdist(XY[train],XY[test]).min(axis=1)>=100]
            assert len(train)>=10 and not set(train)&set(test)
            nearest=cdist(XY[test],XY[train]).min(axis=1)
            foldrecords.append({'scheme':scheme,'fold':int(fold)+1,'train_ids':[int(SITES[i]['id']) for i in train],'test_ids':[int(SITES[i]['id']) for i in test],'n_train':len(train),'n_test':len(test),'minimum_test_to_train_km':float(nearest.min()),'median_test_to_train_km':float(np.median(nearest)),'buffer_km':100 if 'buffer' in scheme else 0})
            configure(train)
            for spec in SPECS:
                try:
                    rows,m=predict(spec,train,test);failure=''
                    fitrecords.append({'scheme':scheme,'fold':int(fold)+1,'model':spec['id'],'boxcox_lambda':m.lam,'parameters':m.params,'primary_trend_coeff':m.primary_coeff.tolist(),'aux_trend_coeff':m.aux_coeff.tolist(),'normalization':m.metadata()['normalization'] if hasattr(m,'metrics') else {'location':m.location,'scale':m.scale},'training_only':True})
                except (ValueError,AssertionError,np.linalg.LinAlgError,RuntimeWarning) as ex:
                    failure=type(ex).__name__+': '+str(ex);rows=[{'station_id':int(SITES[i]['id']),'station':SITES[i]['name'],'observed':float(SITES[i][spec['variable']]),'predicted':None,'se':None,'error':None,'valid':False} for i in test]
                for row in rows:row.update(scheme=scheme,fold=int(fold)+1,model=spec['id'],failure=failure)
                scheme_rows[spec['id']]+=rows
            if scheme!='refit_LOO' or fold%10==0:print('FOLD',scheme,int(fold)+1,'train',len(train),'test',len(test),flush=True)
        for spec in SPECS:
            rows=scheme_rows[spec['id']];assert len(rows)==N and len({r['station_id'] for r in rows})==N
            allmetrics.append(dict(scheme=scheme,model=spec['id'],variable=spec['variable'],method='CK' if spec['co'] else 'OK',**metrics(rows)));allrows+=rows
        csvout(ROOT/'tables'/'validation_metrics.csv',allmetrics);csvout(ROOT/'tables'/'validation_cases.csv',allrows)
    dump(ROOT/'models'/'validation_fits.json',fitrecords);dump(ROOT/'folds.json',foldrecords)
    # Nested 5-fold evaluation: all candidate selection occurs on outer training.
    nestedrows=[];selections=[]
    for outer in np.unique(fold5):
        test=np.flatnonzero(fold5==outer);train=np.flatnonzero(fold5!=outer);inner=groups(XY[train],3,seed=20261101+int(outer))
        innerrows={s['id']:[] for s in SPECS}
        for infold in range(3):
            innertrain=train[inner!=infold];innertest=train[inner==infold];configure(innertrain)
            for spec in SPECS:
                try:rows,m=predict(spec,innertrain,innertest)
                except (ValueError,AssertionError,np.linalg.LinAlgError,RuntimeWarning):rows=[{'valid':False} for i in innertest]
                innerrows[spec['id']]+=rows
        for variable in ['precip','evap']:
            for co in [False,True]:
                candidates=[s for s in SPECS if s['variable']==variable and s['co']==co];scores=[]
                for spec in candidates:
                    rows=innerrows[spec['id']]
                    if len(rows)==len(train) and all(r['valid'] for r in rows):scores.append((float(np.sqrt(np.mean([r['error']**2 for r in rows]))),spec))
                tag=variable+('_CK_nested' if co else '_OK_nested')
                if scores:
                    score,chosen=min(scores,key=lambda r:(r[0],r[1]['id']));configure(train)
                    try:rows,m=predict(chosen,train,test);fail=''
                    except (ValueError,AssertionError,np.linalg.LinAlgError,RuntimeWarning) as ex:rows=[{'station_id':int(SITES[i]['id']),'station':SITES[i]['name'],'observed':float(SITES[i][variable]),'predicted':None,'se':None,'error':None,'valid':False} for i in test];fail=str(ex)
                    selections.append({'outer_fold':int(outer)+1,'variable':variable,'method':'CK' if co else 'OK','chosen':chosen['id'],'inner_RMSE':score,'inner_complete_candidates':len(scores),'outer_n_train':len(train),'outer_n_test':len(test),'training_only_selection':True,'inner_scores':[{'model':s['id'],'RMSE':v} for v,s in sorted(scores,key=lambda x:x[0])]})
                else:raise ValueError('No complete inner candidate: '+tag)
                for r in rows:r.update(scheme='nested_spatial5',fold=int(outer)+1,model=tag,failure=fail)
                nestedrows+=rows
        print('NESTED outer',int(outer)+1,flush=True)
    nestedmetrics=[]
    for variable in ['precip','evap']:
        for co in [False,True]:
            tag=variable+('_CK_nested' if co else '_OK_nested');rows=[r for r in nestedrows if r['model']==tag]
            nestedmetrics.append(dict(scheme='nested_spatial5',model=tag,variable=variable,method='CK' if co else 'OK',**metrics(rows)))
    csvout(ROOT/'tables'/'nested_cases.csv',nestedrows);csvout(ROOT/'tables'/'nested_metrics.csv',nestedmetrics);dump(ROOT/'models'/'nested_selections.json',selections)
    # Parameter identifiability: change range bound AND supported lag bin count.
    sensitivity=[]
    for bound in [1.,2.,4.]:
        for bins in [8,10,14]:
            configure(np.arange(N),{'range_multiplier':bound,'bins':bins})
            for spec in SPECS:
                try:
                    m=model(spec);p=m.params
                    sensitivity.append({'model':spec['id'],'range_multiplier':bound,'bins':bins,'range_km':p['range_km'],'upper_km':p['range_upper_km'],'upper_hit':p['range_upper_hit'],'nugget1':p['nugget1'],'partial1':p['partial1'],'weighted_sse':p['weighted_sse'],'supported_bins':p['supported_bins'],'success':True})
                except Exception as ex:sensitivity.append({'model':spec['id'],'range_multiplier':bound,'bins':bins,'range_km':None,'upper_km':None,'upper_hit':None,'nugget1':None,'partial1':None,'weighted_sse':None,'supported_bins':None,'success':False})
    csvout(ROOT/'tables'/'variogram_sensitivity.csv',sensitivity)
    result={'station_count':N,'fixed_candidate_count':len(SPECS),'metrics':allmetrics,'nested_metrics':nestedmetrics,'nested_selection':'Minimum complete inner spatial3 RMSE, separate OK/CK; outer spatial5 untouched until prediction. Candidate menu was inherited from previous exploratory analysis, not an external holdout.','schemes':list(schemes),'buffer_km':100,'seed':20261030,'folds_depend_only_on_coordinates':True,'all_training_preprocessing_refitted':True,'no_predictions_imputed':True,'constraint_max_error':g.CHECKS['constraint_max_error'],'sensitivity_runs':len(sensitivity)}
    dump(ROOT/'validation_summary.json',result);print('VALIDATION COMPLETE',len(allrows),len(nestedrows),flush=True)

def sensitivity_cv():
    labels=groups(XY,5)
    ids=[PREV[category][variable] for variable in ['precip','evap'] for category in ['ordinary','cokriging']]
    specs=[s for s in SPECS if s['id'] in ids];summaries=[];cases=[]
    for bound in [1.,2.,4.]:
        for bins in [8,10,14]:
            rows={s['id']:[] for s in specs}
            for fold in range(5):
                test=np.flatnonzero(labels==fold);train=np.flatnonzero(labels!=fold);configure(train,{'range_multiplier':bound,'bins':bins})
                for spec in specs:
                    r,m=predict(spec,train,test)
                    for case in r:case.update(model=spec['id'],range_multiplier=bound,bins=bins,fold=fold+1)
                    rows[spec['id']]+=r;cases+=r
            for spec in specs:summaries.append(dict(model=spec['id'],range_multiplier=bound,bins=bins,**metrics(rows[spec['id']])))
            print('Sensitivity spatial5',bound,bins,flush=True)
    csvout(ROOT/'tables'/'sensitivity_spatial5_metrics.csv',summaries);csvout(ROOT/'tables'/'sensitivity_spatial5_cases.csv',cases)
    print('Sensitivity complete',len(summaries))


# Station regression and spatial errors
MAN = {"models": [{"id": dep + "_" + "_".join(terms), "dependent": dep, "predictors": list(terms)}
                  for dep in ("Precip", "Evap") for n_terms in (1, 2, 3)
                  for terms in combinations(("Lat", "Lon", "Elev"), n_terms)]}
FIELD={'Lat':'lat','Lon':'lon','Elev':'elev','Precip':'precip','Evap':'evap'}
def design(spec,ids=None):
    ids=np.arange(N) if ids is None else np.asarray(ids)
    return np.column_stack([np.ones(len(ids))]+[np.array([SITES[i][FIELD[name]] for i in ids]) for name in spec['predictors']])
def yvalues(dep):return np.array([p[FIELD[dep]] for p in SITES])
def moran(e,w):
    e=e-e.mean(axis=0)
    return N/w.sum()*np.sum(e*(w@e),axis=0)/np.sum(e*e,axis=0)
def tails(observed,simulated):
    b=len(simulated);high=float((1+np.sum(simulated>=observed))/(b+1));low=float((1+np.sum(simulated<=observed))/(b+1))
    return high,min(1.,2*min(high,low))
def holm(values):
    order=np.argsort(values);result=np.zeros(len(values));previous=0.
    for rank,i in enumerate(order):previous=max(previous,(len(values)-rank)*values[i]);result[i]=min(1.,previous)
    return result
def fit_sem(x,y,w):
    n=len(y);wx=w@x;wy=w@y;eigen=np.linalg.eigvals(w).real
    def objective(lam,full=False):
        ax=x-lam*wx;ay=y-lam*wy;b=np.linalg.lstsq(ax,ay,rcond=None)[0];e=ay-ax@b;ss=e@e;logdet=np.log(1-lam*eigen).sum();ll=-n/2*(np.log(2*np.pi)+1+np.log(ss/n))+logdet
        return (float(ll),b,e,float(ss/n)) if full else float(-ll)
    grid=np.linspace(-.95,.95,39);values=np.array([objective(a) for a in grid]);i=int(values.argmin());left=grid[max(0,i-1)];right=grid[min(len(grid)-1,i+1)]
    opt=minimize_scalar(objective,bounds=(left,right),method='bounded',options={'xatol':1e-8});lam=float(opt.x);ll,b,e,sigma2=objective(lam,True)
    return {'lambda':lam,'log_likelihood':ll,'beta':b,'innovations':e,'sigma2':sigma2,'boundary_hit':abs(lam)>.949}
def unscale(b,means,scales):return np.r_[b[0]-np.sum(b[1:]*means/scales),b[1:]/scales]
def regression():
    W={}
    dist=cdist(XY,XY);np.fill_diagonal(dist,np.inf)
    for k in [3,4,6,8]:
        a=np.zeros((N,N));a[np.arange(N)[:,None],np.argsort(dist,axis=1)[:,:k]]=1;a=np.maximum(a,a.T);W[k]=a/a.sum(axis=1)[:,None]
    dump(ROOT/'models'/'spatial_weights.json',{'definition':'symmetric union of k-nearest-neighbor edges, then row standardized; projected Euclidean km; no self weights','weights':{str(k):w.tolist() for k,w in W.items()}})
    rng=np.random.default_rng(20261102);B=9999;morans=[];hc3=[];ols={}
    for spec in MAN['models']:
        mid=spec['id'];x=design(spec);y=yvalues(spec['dependent']);q,r=np.linalg.qr(x,mode='reduced');b=np.linalg.lstsq(x,y,rcond=None)[0];e=y-x@b;h=np.sum(q*q,axis=1);df=N-x.shape[1]
        gaussian=rng.normal(size=(N,B));gaussian-=q@(q.T@gaussian)
        wild=e[:,None]/np.sqrt(1-h[:,None])*rng.choice([-1.,1.],size=(N,B));wild-=q@(q.T@wild)
        ri=np.linalg.inv(r);bread=ri@ri.T;cov=bread@(x.T@np.diag((e/(1-h))**2)@x)@bread;se=np.sqrt(np.diag(cov));p=2*stats.t.sf(abs(b/se),df);crit=stats.t.ppf(.975,df)
        for j,name in enumerate(['(Constant)']+spec['predictors']):hc3.append({'model':mid,'term':name,'B':float(b[j]),'HC3_SE':float(se[j]),'HC3_t':float(b[j]/se[j]),'HC3_p_approx':float(p[j]),'CI95_low':float(b[j]-crit*se[j]),'CI95_high':float(b[j]+crit*se[j])})
        ols[mid]={'beta':b.tolist(),'residual':e.tolist()}
        modelrows=[]
        for k,w in W.items():
            observed=float(moran(e,w));sim=moran(gaussian,w);wildsim=moran(wild,w);pg,pt=tails(observed,sim);pw,pwt=tails(observed,wildsim)
            modelrows.append({'model':mid,'neighbors_k':k,'Moran_I':observed,'gaussian_refit_p_positive':pg,'gaussian_refit_p_two_sided':pt,'wild_refit_p_positive':pw,'wild_refit_p_two_sided':pwt,'null_mean_gaussian':float(sim.mean()),'permutations_or_simulations':B})
        adjusted=holm([r['gaussian_refit_p_two_sided'] for r in modelrows])
        for r,p in zip(modelrows,adjusted):r['Holm_p_across_4_weights']=float(p)
        morans+=modelrows
        print('MORAN',mid,modelrows[1]['Moran_I'],modelrows[1]['gaussian_refit_p_two_sided'],modelrows[1]['wild_refit_p_two_sided'],flush=True)
    primary=[r for r in morans if r['neighbors_k']==4 and r['model'] in ['Precip_Lat','Evap_Lat','Precip_Lat_Lon_Elev','Evap_Lat_Lon_Elev']]
    for r,p in zip(primary,holm([r['gaussian_refit_p_two_sided'] for r in primary])):r['primary_Holm_p']=float(p)
    for r in morans:r.setdefault('primary_Holm_p',None)
    csvout(ROOT/'tables'/'residual_moran.csv',morans);csvout(ROOT/'tables'/'HC3_coefficients.csv',hc3)
    # Spatial error model sensitivity, all-data inference only, not claimed CV.
    sem=[];bootstrap_records=[]
    for dep in ['Precip','Evap']:
        spec=next(s for s in MAN['models'] if s['id']==dep+'_Lat_Lon_Elev');xraw=design(spec);means=xraw[:,1:].mean(axis=0);scales=xraw[:,1:].std(axis=0,ddof=1);x=np.column_stack([np.ones(N),(xraw[:,1:]-means)/scales]);y=yvalues(dep)
        for k in [4,8]:
            w=W[k];fit=fit_sem(x,y,w);a=np.eye(N)-fit['lambda']*w;ainv=np.linalg.inv(a);innovation=fit['innovations']-fit['innovations'].mean();center=x@fit['beta'];rep=[];rep_lambda=[];rep_moran=[];boundaries=0
            for boot in range(499):
                # Empirical centered innovation bootstrap preserves non-Gaussian shape.
                sampled=rng.choice(innovation,size=N,replace=True);yb=center+ainv@sampled;fb=fit_sem(x,yb,w);rep.append(unscale(fb['beta'],means,scales));rep_lambda.append(fb['lambda']);rep_moran.append(float(moran(fb['innovations'],w)));boundaries+=int(fb['boundary_hit'])
            rep=np.array(rep);beta=unscale(fit['beta'],means,scales);rawres=y-x@fit['beta'];observedI=float(moran(fit['innovations'],w));_,ip=tails(observedI,np.array(rep_moran));olsfit=fit_sem(x,y,np.zeros_like(w));lr=max(0.,2*(fit['log_likelihood']-olsfit['log_likelihood']))
            coefficients=[]
            for j,name in enumerate(['(Constant)']+spec['predictors']):
                low,high=np.quantile(rep[:,j],[.025,.975]);coefficients.append({'term':name,'B':float(beta[j]),'bootstrap_SE':float(rep[:,j].std(ddof=1)),'percentile_CI95_low':float(low),'percentile_CI95_high':float(high)})
            row={'model':spec['id'],'neighbors_k':k,'lambda':fit['lambda'],'lambda_bootstrap_CI95':np.quantile(rep_lambda,[.025,.975]).tolist(),'lambda_boundary_hit':fit['boundary_hit'],'log_likelihood':fit['log_likelihood'],'LR_vs_OLS':lr,'LR_p_chi_square_approx':float(stats.chi2.sf(lr,1)),'innovation_Moran_I':observedI,'innovation_refitted_bootstrap_p_two_sided':ip,'bootstrap_replicates':499,'bootstrap_boundary_fraction':boundaries/499,'coefficients':coefficients,'raw_residual':rawres.tolist(),'innovations':fit['innovations'].tolist(),'note':'SEM y=Xb+u, u=lambda Wu+epsilon. Empirical-innovation bootstrap is conditional on this spatial error specification; no claim of spatially valid prediction accuracy.'}
            sem.append(row);bootstrap_records.append({'model':spec['id'],'k':k,'lambda':rep_lambda,'beta':rep.tolist(),'innovation_I':rep_moran})
            print('SEM',dep,'k',k,'lambda',fit['lambda'],'CI',row['lambda_bootstrap_CI95'],'I_p',ip,flush=True)
    dump(ROOT/'models'/'SEM_models.json',sem);dump(ROOT/'models'/'SEM_bootstrap.json',bootstrap_records)
    # OLS shares exactly the geographic folds/buffer used for kriging.
    folds5=groups(XY,5);folds8=groups(XY,8);schemes={'refit_LOO':np.arange(N),'spatial5':folds5,'spatial8':folds8,'spatial5_buffer100km':folds5};cvrows=[];cvmetrics=[]
    for scheme,labels in schemes.items():
        for spec in MAN['models']:
            x=design(spec);y=yvalues(spec['dependent']);rows=[]
            for fold in np.unique(labels):
                test=np.flatnonzero(labels==fold);train=np.flatnonzero(labels!=fold)
                if 'buffer' in scheme:train=train[cdist(XY[train],XY[test]).min(axis=1)>=100]
                beta=np.linalg.lstsq(x[train],y[train],rcond=None)[0];pred=x[test]@beta
                for i,p in zip(test,pred):rows.append({'scheme':scheme,'model':spec['id'],'fold':int(fold)+1,'station_id':int(SITES[i]['id']),'observed':float(y[i]),'predicted':float(p),'error':float(p-y[i])})
            error=np.array([r['error'] for r in rows]);cvmetrics.append({'scheme':scheme,'model':spec['id'],'variable':FIELD[spec['dependent']],'n_valid':len(rows),'RMSE':float(np.sqrt(np.mean(error**2))),'ME':float(error.mean())});cvrows+=rows
    csvout(ROOT/'tables'/'OLS_spatial_cv.csv',cvmetrics);csvout(ROOT/'tables'/'OLS_spatial_cases.csv',cvrows)
    nested=[];choices=[]
    for outer in range(5):
        test=np.flatnonzero(folds5==outer);train=np.flatnonzero(folds5!=outer);inner=groups(XY[train],3,seed=20261101+outer)
        for dep in ['Precip','Evap']:
            scores=[]
            for spec in [s for s in MAN['models'] if s['dependent']==dep]:
                x=design(spec);y=yvalues(dep);errors=[]
                for fold in range(3):
                    itrain=train[inner!=fold];itest=train[inner==fold];beta=np.linalg.lstsq(x[itrain],y[itrain],rcond=None)[0];errors.extend((x[itest]@beta-y[itest]).tolist())
                scores.append((float(np.sqrt(np.mean(np.array(errors)**2))),spec))
            score,chosen=min(scores,key=lambda p:(p[0],p[1]['id']));x=design(chosen);y=yvalues(dep);b=np.linalg.lstsq(x[train],y[train],rcond=None)[0];pred=x[test]@b
            choices.append({'outer_fold':outer+1,'dependent':dep,'chosen':chosen['id'],'inner_RMSE':score})
            for i,p in zip(test,pred):nested.append({'variable':FIELD[dep],'model':FIELD[dep]+'_OLS_nested','station_id':int(SITES[i]['id']),'fold':outer+1,'observed':float(y[i]),'predicted':float(p),'error':float(p-y[i])})
    nestedmetrics=[]
    for variable in ['precip','evap']:
        rows=[r for r in nested if r['variable']==variable];err=np.array([r['error'] for r in rows]);nestedmetrics.append({'scheme':'nested_spatial5','variable':variable,'model':variable+'_OLS_nested','n_valid':len(rows),'RMSE':float(np.sqrt(np.mean(err**2))),'ME':float(err.mean())})
    csvout(ROOT/'tables'/'OLS_nested_cases.csv',nested);dump(ROOT/'models'/'OLS_nested_choices.json',choices)
    dump(ROOT/'regression_spatial_summary.json',{'n':N,'Moran_null':'9999 simulated iid Gaussian innovations projected through residual-maker I-H, plus heteroskedasticity-sensitive HC2-scaled residual wild bootstrap refit; two-sided tail p. Not naive permutations of constrained OLS residuals.','weights':'Symmetric union kNN then row standardization; k=3,4,6,8; k4 primary; no islands','primary_Moran':primary,'all_Moran':morans,'HC3_note':'HC3 allows heteroskedasticity, not spatial correlation. Approximate t intervals.','SEM':sem,'OLS_spatial_metrics':cvmetrics,'OLS_nested_metrics':nestedmetrics,'SEM_note':'Spatial-error correction is a sensitivity analysis of full linear-model inference, not a replacement for all geography or proof of independent innovations.'})

if __name__ == "__main__":
    main()
    sensitivity_cv()
    regression()
