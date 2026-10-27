# print.galick_run snapshot for successful run

    Code
      print(r)
    Message
      
      -- galick_run: OK --------------------------------------------------------------
      Runtime: 3.45 sec
      Command: java -jar beagle.jar gt=in.vcf out=out
      
      -- Outputs 
      * '/tmp/out.vcf.gz'
      
      -- Logs 
      * stdout: '/tmp/stdout.log'
      * stderr: '/tmp/stderr.log'

# print.galick_run snapshot for failed run

    Code
      print(r)
    Message
      
      -- galick_run: status 7 --------------------------------------------------------
      Runtime: 0.50 sec
      Command: java -jar beagle.jar
      
      -- Logs 
      * stderr: '/tmp/stderr.log'

