# -*- coding: utf-8 -*-
"""Ordinary kriging and genuine two-variable ordinary cokriging (LMC).

Guide 5(1): data/preprocessing -> detrending -> variogram -> neighbors -> CV
-> model summary -> prediction. Guide 5(2): primary + station elevation.
The analysis grid uses projected cell centres inside Gansu.
"""
import os
for env in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[env]='1'
import json,csv,warnings
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist,pdist,squareform
from scipy.linalg import solve,LinAlgWarning
from scipy.special import boxcox,inv_boxcox
import shapely
import rasterio
from rasterio.transform import from_origin

from scipy import stats
from inputs import WORK, station_inputs

ROOT = WORK
for folder in ["models", "tables", "rasters"]:
    (ROOT / folder).mkdir(parents=True, exist_ok=True)
SITES, BOUNDARY = station_inputs()
N = len(SITES)
XY = np.array([p["xy"] for p in SITES]) / 1000
ESDA = {"variables": {v: {"boxcox_lambda_mle": float(stats.boxcox([p[v] for p in SITES])[1])}
                      for v in ("precip", "evap")}}
ORIGIN=XY.mean(axis=0);AXIS_SCALE=XY.std(axis=0,ddof=1)
II,JJ=np.triu_indices(N,1);PAIR_DIST=pdist(XY);D=squareform(PAIR_DIST)
MAXLAG=float(PAIR_DIST.max()/2);EDGES=np.linspace(0,MAXLAG,11)
if not BOUNDARY.is_valid:BOUNDARY=shapely.make_valid(BOUNDARY)
BBOX=BOUNDARY.bounds
RADIUS=float(np.hypot(max(BBOX[2]/1000,XY[:,0].max())-min(BBOX[0]/1000,XY[:,0].min()),max(BBOX[3]/1000,XY[:,1].max())-min(BBOX[1]/1000,XY[:,1].min())))
CELL=5000.;XMIN=np.floor(BBOX[0]/CELL)*CELL;XMAX=np.ceil(BBOX[2]/CELL)*CELL
YMIN=np.floor(BBOX[1]/CELL)*CELL;YMAX=np.ceil(BBOX[3]/CELL)*CELL
GX=np.arange(XMIN+CELL/2,XMAX,CELL);GY=np.arange(YMAX-CELL/2,YMIN,-CELL)
XX,YY=np.meshgrid(GX,GY);MASK=shapely.contains_xy(BOUNDARY,XX,YY)
GRID=np.column_stack([XX[MASK],YY[MASK]])/1000
FIT_CACHE={};CHECKS={'n':N,'pair_count':len(PAIR_DIST),'constraint_max_error':0.,'covariance_min_eigenvalues':{},'jitter_normalized':1e-10}

def write_json(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def write_csv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,list(rows[0]));writer.writeheader();writer.writerows(rows)
def design(q):
    x,y=((np.asarray(q)-ORIGIN)/AXIS_SCALE).T
    return np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
def corr(h,kind,a):
    t=np.asarray(h)/a
    if kind=='spherical':return np.where(t<1,1-1.5*t+.5*t**3,0.)
    if kind=='exponential':return np.exp(-3*t)
    if kind=='gaussian':return np.exp(-3*t*t)
    raise ValueError(kind)
def empirical(u,a=None):
    du=u[II]-u[JJ];da=None if a is None else a[II]-a[JJ]
    rows=[]
    for k in range(10):
        use=(PAIR_DIST>=EDGES[k])&(PAIR_DIST<=EDGES[k+1] if k==9 else PAIR_DIST<EDGES[k+1])
        if use.sum()<5:continue
        rows.append(dict(h=float(PAIR_DIST[use].mean()),count=int(use.sum()),g1=float((.5*du[use]**2).mean()),g2=None if a is None else float((.5*da[use]**2).mean()),g12=None if a is None else float((.5*du[use]*da[use]).mean())))
    assert len(rows)>=5
    return rows
def fit_variogram(u,kind,aux=None):
    emp=empirical(u,aux);h=np.array([r['h'] for r in emp]);weight=np.sqrt([r['count'] for r in emp]);weight/=weight.max()
    g1=np.array([r['g1'] for r in emp]);g2=None if aux is None else np.array([r['g2'] for r in emp]);g12=None if aux is None else np.array([r['g12'] for r in emp])
    minrange=max(5.,float(np.min(PAIR_DIST)/2));maxrange=float(PAIR_DIST.max()*2)
    if aux is None:
        lo=[minrange,1e-6,1e-8];hi=[maxrange,10.,3.]
        def residual(p):return weight*(p[2]+p[1]*(1-corr(h,kind,p[0]))-g1)
    else:
        lo=[minrange,1e-6,1e-6,1e-8,1e-8,-.995];hi=[maxrange,10.,10.,3.,3.,.995]
        def residual(p):
            shape=1-corr(h,kind,p[0]);cross=p[5]*np.sqrt(p[1]*p[2])
            return np.concatenate([weight*(p[3]+p[1]*shape-g1),weight*(p[4]+p[2]*shape-g2),weight*(cross*shape-g12)])
    fits=[]
    for factor in [.2,.6,1.2]:
        r0=float(np.clip(PAIR_DIST.max()*factor,minrange*1.01,maxrange*.99))
        x0=[r0,1.,.05] if aux is None else [r0,1.,1.,.05,.05,float(np.clip(np.corrcoef(u,aux)[0,1],-.8,.8))]
        fit=least_squares(residual,x0,bounds=(lo,hi),loss='linear',max_nfev=3000,x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10)
        if fit.success and np.isfinite(fit.x).all():fits.append(fit)
    assert fits,'Variogram fitting did not converge'
    best=min(fits,key=lambda r:np.sum(r.fun**2));p=best.x
    result={'range_km':float(p[0]),'partial1':float(p[1]),'nugget1':float(p[2] if aux is None else p[3]),'weighted_sse':float(np.sum(best.fun**2)),'bound_hits':[i for i,v in enumerate(p) if abs(v-lo[i])<1e-5 or abs(v-hi[i])<1e-5]}
    if aux is not None:result.update(partial2=float(p[2]),nugget2=float(p[4]),rho=float(p[5]),cross_partial=float(p[5]*np.sqrt(p[1]*p[2])))
    return result,emp

class Model:
    def __init__(self,variable,kind='spherical',transform='raw',trend=0,k=12,co=False):
        self.variable=variable;self.kind=kind;self.transform=transform;self.trend=trend;self.k=int(k);self.co=co
        self.id='%s_%s_%s_%s_T%d_N%d'%(variable,'CK' if co else 'OK',transform,kind,trend,k)
        self.z=np.array([r[variable] for r in SITES],float)
        self.lam=float(ESDA['variables'][variable]['boxcox_lambda_mle']) if transform=='boxcox' else None
        value=self.z if transform=='raw' else np.log(self.z) if transform=='log' else boxcox(self.z,self.lam)
        self.location=float(value.mean());self.scale=float(value.std(ddof=1));self.u=(value-self.location)/self.scale
        elevation=np.array([r['elev'] for r in SITES],float);self.aux_location=float(elevation.mean());self.aux_scale=float(elevation.std(ddof=1));self.a=(elevation-self.aux_location)/self.aux_scale
        self.primary_coeff=np.linalg.lstsq(design(XY),self.u,rcond=None)[0] if trend==2 else np.zeros(6)
        self.aux_coeff=np.linalg.lstsq(design(XY),self.a,rcond=None)[0] if trend==2 else np.zeros(6)
        residual=self.u-design(XY)@self.primary_coeff;aux=self.a-design(XY)@self.aux_coeff
        cache=(variable,kind,transform,trend,co)
        if cache not in FIT_CACHE:FIT_CACHE[cache]=fit_variogram(residual,kind,aux if co else None)
        self.params,self.empirical=FIT_CACHE[cache]
        c=self.covariance(D)
        eig=float(np.linalg.eigvalsh(c).min());CHECKS['covariance_min_eigenvalues'][self.id]=eig
        assert eig>-1e-8,'Invalid joint covariance'
    def covariance(self,h):
        p=self.params;r=corr(h,self.kind,p['range_km']);zero=np.asarray(h)<1e-9
        c11=p['partial1']*r+p['nugget1']*zero
        if not self.co:return c11
        c22=p['partial2']*r+p['nugget2']*zero;c12=p['cross_partial']*r
        return np.block([[c11,c12],[c12,c22]])
    def inverse(self,mu,sigma):
        y=mu*self.scale+self.location;sd=sigma*self.scale
        with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
            if self.transform=='raw':return y,sd
            if self.transform=='log':out=np.exp(y);return out,out*sd
            out=inv_boxcox(y,self.lam)
            derivative=np.power(self.lam*y+1,1/self.lam-1) if abs(self.lam)>1e-8 else np.exp(y)
            return out,np.abs(derivative)*sd
    def predict(self,q,exclude=None):
        q=np.atleast_2d(q);available=np.array([i for i in range(N) if i!=exclude]);k=min(self.k,len(available))
        if exclude is not None:assert exclude not in available
        dist=cdist(q,XY[available]);nearest=available[np.argsort(dist,axis=1,kind='stable')[:,:k]]
        nearest.sort(axis=1);groups,inverse=np.unique(nearest,axis=0,return_inverse=True)
        if self.trend==2 and exclude is not None:
            cp=np.linalg.lstsq(design(XY[available]),self.u[available],rcond=None)[0];ca=np.linalg.lstsq(design(XY[available]),self.a[available],rcond=None)[0]
        else:cp=self.primary_coeff;ca=self.aux_coeff
        u=self.u-design(XY)@cp;a=self.a-design(XY)@ca
        pred=np.empty(len(q));variance=np.empty(len(q))
        for gi,idx in enumerate(groups):
            at=np.flatnonzero(inverse==gi);h=cdist(XY[idx],q[at]);assert h.max()<=RADIUS+1e-6
            p=self.params;core=corr(h,self.kind,p['range_km'])
            cov=self.covariance(D[np.ix_(idx,idx)]);cov=cov+np.eye(len(cov))*1e-10
            if self.co:
                f=np.zeros((2*k,2));f[:k,0]=1;f[k:,1]=1
                rhs_cov=np.vstack([p['partial1']*core+p['nugget1']*(h<1e-9),p['cross_partial']*core])
                rhs=np.vstack([rhs_cov,np.ones((1,len(at))),np.zeros((1,len(at)))]);obs=np.concatenate([u[idx],a[idx]])
                target_constraint=np.array([1.,0.])[:,None]
            else:
                f=np.ones((k,1));rhs_cov=p['partial1']*core+p['nugget1']*(h<1e-9)
                rhs=np.vstack([rhs_cov,np.ones((1,len(at)))]);obs=u[idx];target_constraint=np.ones((1,1))
            aug=np.block([[cov,f],[f.T,np.zeros((f.shape[1],f.shape[1]))]])
            with warnings.catch_warnings():
                warnings.simplefilter('error',LinAlgWarning)
                weights=solve(aug,rhs,assume_a='sym')
            w=weights[:len(obs)];mu=weights[len(obs):]
            constraint=float(np.max(np.abs(f.T@w-target_constraint)));CHECKS['constraint_max_error']=max(CHECKS['constraint_max_error'],constraint)
            assert constraint<1e-6,'Unbiasedness constraints failed'
            pred[at]=obs@w+design(q[at])@cp
            variance[at]=p['partial1']+p['nugget1']-np.sum(w*rhs_cov,axis=0)-mu[0]
        assert variance.min()>-1e-6,'Negative kriging variance'
        return self.inverse(pred,np.sqrt(np.maximum(variance,0.)))
    def cv(self):
        rows=[]
        for i in range(N):
            pred,se=self.predict(XY[i],exclude=i);p=float(pred[0]);s=float(se[0]);valid=bool(np.isfinite(p) and np.isfinite(s) and s>0)
            rows.append(dict(model=self.id,station_id=int(SITES[i]['id']),station=SITES[i]['name'],observed=float(self.z[i]),predicted=p if valid else '',standard_error=s if valid else '',error=p-self.z[i] if valid else '',standardized_error=(p-self.z[i])/s if valid else '',valid=valid))
        valid=[r for r in rows if r['valid']]
        metrics={'model':self.id,'variable':self.variable,'method':'CK' if self.co else 'OK','transform':self.transform,'trend':self.trend,'variogram':self.kind,'neighbors':self.k,'n_valid':len(valid)}
        if len(valid)==N:
            err=np.array([r['error'] for r in rows]);se=np.array([r['standard_error'] for r in rows]);zerr=err/se
            metrics.update(ME=float(err.mean()),RMSE=float(np.sqrt(np.mean(err**2))),ASE=float(se.mean()),MeanStandardized=float(zerr.mean()),RMSSE=float(np.sqrt(np.mean(zerr**2))))
        else:metrics.update(ME=None,RMSE=None,ASE=None,MeanStandardized=None,RMSSE=None)
        self.cv_rows=rows;self.metrics=metrics
        write_csv(ROOT/'tables'/(self.id+'_cv.csv'),rows)
        return metrics
    def metadata(self):
        return dict(id=self.id,variable=self.variable,method='ordinary cokriging' if self.co else 'ordinary kriging',auxiliary='station elevation' if self.co else None,transform=self.transform,boxcox_lambda=self.lam,trend_order=self.trend,neighbors=self.k,min_neighbors=3,sector_count=1,search_radius_km=RADIUS,kernel=self.kind,parameters=self.params,empirical_bins=self.empirical,normalization={'location':self.location,'scale':self.scale,'aux_location':self.aux_location,'aux_scale':self.aux_scale},primary_trend_coeff=self.primary_coeff.tolist(),aux_trend_coeff=self.aux_coeff.tolist(),xy_basis_origin_km=ORIGIN.tolist(),xy_basis_scale_km=AXIS_SCALE.tolist(),crossvalidation=self.metrics,backtransform='inverse-transform central prediction (median on original scale); standard error via first-order delta method',cv_support='Leave one station out of BOTH primary and elevation. Variogram and transform parameters fitted on full dataset; polynomial trend refitted on each training fold. Conditional model diagnostics, not independent external validation.')

def select(models):
    valid=[m for m in models if m.metrics['n_valid']==N and m.metrics['RMSE'] is not None and m.metrics.get('grid_missing',0)==0]
    assert valid,'No complete candidate'
    best=min(m.metrics['RMSE'] for m in valid)
    close=[m for m in valid if m.metrics['RMSE']<=best*1.01]
    return min(close,key=lambda m:(abs(m.metrics['RMSSE']-1),m.trend,m.transform!='raw',m.metrics['RMSE'],m.id))
def map_model(model):
    pred,se=model.predict(GRID);good=np.isfinite(pred)&np.isfinite(se)&(se>=0)
    model.metrics['grid_missing']=int((~good).sum());model.metrics['grid_negative']=int(np.sum(pred[good]<0));model.metrics['grid_cells']=int(MASK.sum())
    for suffix,value in [('prediction',pred),('se',se)]:
        array=np.full(MASK.shape,-9999,dtype=np.float32);arr=np.where(good,value,-9999).astype(np.float32);array[MASK]=arr
        filename=ROOT/'rasters'/(model.id+'_'+suffix+'.tif')
        with rasterio.open(filename,'w',driver='GTiff',height=array.shape[0],width=array.shape[1],count=1,dtype='float32',crs='EPSG:2381',transform=from_origin(XMIN,YMAX,CELL,CELL),nodata=-9999,compress='deflate') as dst:
            dst.write(array,1);dst.update_tags(model_id=model.id,algorithm='script ordinary cokriging LMC' if model.co else 'script ordinary kriging',standard_error='conditional; delta approximation after nonlinear transform',pixel_interpretation='point prediction at cell center')
        with rasterio.open(filename) as src:
            assert src.crs.to_epsg()==2381 and src.res==(5000.,5000.)
            reread=src.read(1);assert np.array_equal(reread!=-9999,MASK&np.isfinite(array)&(array!=-9999))
    write_json(ROOT/'models'/(model.id+'.json'),model.metadata())
    print('GRID',model.id,'missing',model.metrics['grid_missing'],'negative',model.metrics['grid_negative'],flush=True)
def verify_solver():
    h=np.array([[0.,2.],[2.,0.]]);c=corr(h,'spherical',3.)
    rhs=np.r_[corr(np.array([1.,1.]),'spherical',3.),1.]
    weights=np.linalg.solve(np.block([[c,np.ones((2,1))],[np.ones((1,2)),np.zeros((1,1))]]),rhs)
    assert np.allclose(weights[:2],[.5,.5],atol=1e-12)
    assert np.isclose(weights[:2]@np.array([2.,4.]),3.)
    # Independent secondary variable must have zero influence and reduce to OK.
    full=np.block([[c,np.zeros((2,2))],[np.zeros((2,2)),c]])
    f=np.array([[1,0],[1,0],[0,1],[0,1]],float)
    co=np.linalg.solve(np.block([[full,f],[f.T,np.zeros((2,2))]]),np.r_[rhs[:2],[0.,0.],[1.,0.]])
    assert np.allclose(co[:2],weights[:2],atol=1e-12) and np.allclose(co[2:4],0,atol=1e-12)
    CHECKS['symmetric_two_point_prediction']=3.;CHECKS['zero_crosscovariance_reduces_to_OK']=True

def main():
    verify_solver();all_models=[];best_ok={};baselines={};best_ck={}
    def evaluate(model):
        model.cv();print('CV',model.id,'RMSE',model.metrics['RMSE'],'RMSSE',model.metrics['RMSSE'],flush=True)
        all_models.append(model);write_json(ROOT/'models'/(model.id+'.json'),model.metadata());return model
    # Guide 5(1): finish both single-variable workflows before any cokriging.
    for variable in ['precip','evap']:
        family=[evaluate(Model(variable,kind=kind)) for kind in ['spherical','exponential','gaussian']]
        baselines[variable]=family[0];base=select(family)
        transforms=[base]+[evaluate(Model(variable,kind=base.kind,transform=t)) for t in ['log','boxcox']]
        chosen=select(transforms)
        trend=evaluate(Model(variable,kind=chosen.kind,transform=chosen.transform,trend=2))
        chosen=select([chosen,trend])
        for k in [8,16]:evaluate(Model(variable,kind=chosen.kind,transform=chosen.transform,trend=chosen.trend,k=k))
        candidates=[m for m in all_models if m.variable==variable and not m.co]
        for model in candidates:map_model(model)
        best_ok[variable]=select(candidates)
        print('BEST_OK',variable,best_ok[variable].id,flush=True)
    # Guide 5(2): actual bivariate ordinary cokriging, not elevation regression.
    for variable in ['precip','evap']:
        base=best_ok[variable]
        choices = []
        for kind in ("spherical", "exponential", "gaussian"):
            model = evaluate(Model(variable, kind=kind, transform=base.transform,
                                   trend=base.trend, k=base.k, co=True))
            map_model(model)
            choices.append(model)
        best_ck[variable] = select(choices)
    # Exact-observation check on the untransformed default, numerical jitter only.
    for variable,m in baselines.items():
        pred,se=m.predict(XY)
        assert np.allclose(pred,m.z,rtol=1e-6,atol=1e-4)
        CHECKS[variable+'_observed_value_reproduction_max_error']=float(np.max(np.abs(pred-m.z)))
    for m in all_models:write_json(ROOT/'models'/(m.id+'.json'),m.metadata())
    write_csv(ROOT/'tables'/'model_comparison.csv',[m.metrics for m in all_models])
    write_json(ROOT/'solver_checks.json',CHECKS)
    summary={'engine':'NumPy/SciPy ordinary kriging and LMC ordinary cokriging; not claimed as ArcGIS Geostatistical Analyst output','n':N,'models':len(all_models),'selection_rule':'53 complete CV sites; minimum RMSE, within 1% prefer RMSSE closer to 1, then simpler preprocessing; exclude incomplete raster domain','grid':{'cell_size_m':CELL,'width':len(GX),'height':len(GY),'inside_cells':int(MASK.sum()),'xmin':XMIN,'ymin':YMIN,'xmax':XMAX,'ymax':YMAX,'epsg':2381},'variogram_bins':{'count':10,'maxlag_km':MAXLAG,'weight':'square root of pair count; weighted least squares; same bins for all models'},'baseline_definition':'Script baseline: raw values, spherical variogram, no detrending, 12 nearest sites, one sector. Not ArcGIS wizard defaults.','cv_method':'Leave-one-STATION-out (primary and elevation both removed); variogram and transformation parameters fitted once using full data; polynomial refitted per fold. Error = predicted - observed.','uncertainty_note':'Conditional on fitted variogram/trend; nonlinear backtransform uses inverse central prediction and delta-method standard error. Parameter/trend-estimation uncertainty not fully included.','cokriging_model':'LMC: structured PSD B1 with rho bounded [-.995,.995], diagonal nonnegative nugget B0; separate ordinary-mean constraints primary sum=1, elevation sum=0. Only station elevation observations, no DEM assumed.','ordinary':{v:best_ok[v].id for v in best_ok},'cokriging':{v:best_ck[v].id for v in best_ck},'baseline':{v:baselines[v].id for v in baselines},'model_metrics':[m.metrics for m in all_models]}
    write_json(ROOT/'summary.json',summary)
    print(json.dumps({'models':len(all_models),'best_ordinary':summary['ordinary'],'cokriging':summary['cokriging'],'checks':{k:v for k,v in CHECKS.items() if k!='covariance_min_eigenvalues'}},ensure_ascii=True,indent=2),flush=True)

if __name__=='__main__':main()
