from .liljegren import liljegren as liljegrenWBGT
from .bernard import bernard as bernardWBGT
from .dimiceli import dimiceli as dimiceliWBGT

def wbgt( method, *args, **kwargs ):

  method = method.lower()
  if method == 'liljegren':
    return liljegrenWBGT( *args, **kwargs )
  elif method == 'bernard':
    return bernardWBGT( *args, **kwargs )
  elif method == 'dimiceli':
    return dimiceliWBGT( *args, **kwargs )

  raise Exception( f'Unsupported WBGT method : {method}' )
     
