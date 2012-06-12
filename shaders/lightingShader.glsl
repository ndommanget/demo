// lightingShader.glsl
// Fluids Simulation
// N. Dommanget ndommanget@gmail.com



#ifdef _VERTEX_

// Attributes : per vertex data
in vec4 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexUvs;
in vec4 vertexColor;

// Varyings : data to transmit to fragments
smooth out vec4 position;
smooth out vec4 normal;
smooth out vec2 uvs;
smooth out vec4 localColor;

void main()
{
    if (filledData[0]) position = model * vertexPosition;
    if (filledData[1]) normal = model * vec4(vertexNormal, 0.0);
    if (filledData[2]) uvs = vertexUvs;
    if (filledData[3]) localColor = vertexColor;

    gl_Position = projection * view * model * vertexPosition;
}


#endif




#ifdef _FRAGMENT_


// Varyings : data receved and interpolated from the vertex shaders
smooth in vec4 position;
smooth in vec4 normal;
smooth in vec2 uvs;
smooth in vec4 localColor;

// Final output
out vec4 fragColor;

void main()
{
    vec4 diffuse=vec4(material.diffuse);
    // If color
    if (filledData[3]) diffuse=localColor;
    
    // If no normal
    if (!filledData[1]) fragColor=vec4(diffuse.rgb*diffuse.a, diffuse.a);
    else
    {
        vec4 L=normalize(light.position); // Direction of light from fragment -> light.position[3]==0.0 : Directional light
        if (light.position.w==1.0) L=normalize(light.position-position); //   -> light.position[3]==1.0 : Point light
        vec4 V=normalize(center-position); // Direction from fragment to camera center
        vec4 R=normalize(reflect(-L, normal)); // Direction of reflected light beam, from fragment
        vec4 N=normalize(normal); // Normal
    
        float lightPower=light.power;
        // If spot to the center
        vec4 spotLightDir=vec4(normalize(vec3(0.0, 0.0, 0.0)-light.position.xyz), 1.0);
        float dot=dot(spotLightDir, -L);
        float cutOffMax=0.99;
        float cutOffMin=0.92;
        if (dot<cutOffMax)
        { 
            float decay=0.3;
            float coef=decay;
            if (dot>cutOffMin)
                coef=(0.5*cos((cutOffMax-dot)/(cutOffMax-cutOffMin)*3.14)-0.5)*(1.0-decay)+1.0;
            lightPower=coef*lightPower;
        }

        float ambientValue=lightPower;
        float diffuseValue=lightPower * max( dot(L, N), 0.0);
        float specularValue=lightPower * pow( max(dot(R, V), 0.0), material.shininess);

        vec3 ambientContribution=ambientValue*material.ka*material.ambient.rgb*material.ambient.a;
        vec3 diffuseContribution=diffuseValue*material.kd*diffuse.rgb*diffuse.a;
        vec3 specularContribution=specularValue*material.ks*material.specular.rgb*material.specular.a;

        float alpha=(material.ambient.a*material.ka+diffuse.a*material.kd+material.specular.a*material.ks)/(material.ka+material.kd+material.ks);

        fragColor=vec4(ambientContribution+diffuseContribution+specularContribution, alpha);
    }
}


#endif
